

rule uk_strip_header_digits:
    input:
        fasta = config["latest_majora_matched_fasta"],
    output:
        fasta = config["output_path"] + "/1a/uk_latest.headerstripped.fasta",
    log:
        config["output_path"] + "/logs/1a_uk_strip_header_digits.log"
    run:
        from Bio import SeqIO

        fasta_in = SeqIO.parse(str(input.fasta), "fasta")
        with open(str(output.fasta), 'w') as f:
            for record in fasta_in:
                ID = record.description.split("|")[0]
                f.write(">" + ID + "\n")
                f.write(str(record.seq) + "\n")


rule uk_add_sample_date:
    input:
        metadata = config["latest_majora_matched_metadata"],
    output:
        metadata = config["output_path"] + "/1a/uk_latest.add_sample_date.csv",
    log:
        config["output_path"] + "/logs/1a_uk_add_sample_date.log"
    run:
        from datetime import datetime
        import csv

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
        open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter="\t", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["sample_date"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                try:
                    date = datetime.strptime(row["collection_date"], '%Y-%m-%d').date()
                    row["sample_date"] = row["collection_date"]
                except:
                    try:
                        date = datetime.strptime(row["received_date"], '%Y-%m-%d').date()
                        row["sample_date"] = row["received_date"]
                    except:
                        row["sample_date"] = ""

                writer.writerow(row)


rule uk_update_sample_dates:
    input:
        metadata = rules.uk_add_sample_date.output.metadata,
        updated_dates = config["uk_updated_dates"],
    output:
        metadata = config["output_path"] + "/1a/uk_latest.update_sample_date.csv",
    log:
        config["output_path"] + "/logs/1a_uk_update_sample_dates.log"
    resources: mem_per_cpu=5000
    run:
        import csv

        date_dict = {}

        with open(str(input.updated_dates), 'r', newline = '') as dates_in:
            reader = csv.DictReader(dates_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:
                date_dict[row["central_sample_id"]] = row["sample_date"]


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                if row["central_sample_id"] in date_dict:
                    row["sample_date"] = date_dict[row["central_sample_id"]]

                writer.writerow(row)


rule uk_filter_on_sample_date:
    input:
        metadata = rules.uk_update_sample_dates.output.metadata,
        fasta = rules.uk_strip_header_digits.output.fasta,
    params:
        time_window = config["time_window"],
        date = config["date"],
    output:
        metadata = config["output_path"] + "/1a/uk_latest.filter_on_days.csv",
        fasta = config["output_path"] + "/1a/uk_latest.filter_on_days.fasta",
    log:
        config["output_path"] + "/logs/1a_uk_filter_on_sample_date.log"
    run:
        import datetime
        from Bio import SeqIO
        import csv

        indexed_fasta = SeqIO.index(str(input.fasta), "fasta")

        window = datetime.timedelta(int(params.time_window))
        todays_date = datetime.datetime.strptime(str(params.date), '%Y-%m-%d').date()

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.fasta), "w") as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                try:
                    date = datetime.datetime.strptime(row["sample_date"], '%Y-%m-%d').date()
                except:
                    continue

                if (todays_date - window) > date:
                    continue

                writer.writerow(row)

                seq_rec = indexed_fasta[row["fasta_header"]]
                fasta_out.write(">" + seq_rec.id + "\n")
                fasta_out.write(str(seq_rec.seq) + "\n")


rule uk_filter_omitted_sequences:
    input:
        fasta = rules.uk_filter_on_sample_date.output.fasta,
        metadata = rules.uk_filter_on_sample_date.output.metadata,
        omissions = config["uk_omissions"]
    output:
        fasta = config["output_path"] + "/1a/uk_latest.filter_on_omitted.fasta",
        metadata = config["output_path"] + "/1a/uk_latest.filter_on_omitted.csv",
    log:
        config["output_path"] + "/logs/1a_uk_filter_omitted_sequences.log"
    run:
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index(str(input.fasta), "fasta")

        omissions = set()

        with open(str(input.omissions), "r") as f:
            for line in f:
                omissions.add(line.rstrip())

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.fasta), "w") as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                if row["central_sample_id"] in omissions:
                    continue

                record = alignment[row["fasta_header"]]

                writer.writerow(row)
                fasta_out.write(">" + record.id + "\n")
                fasta_out.write(str(record.seq) + "\n")


rule uk_add_previous_pangolin_lineages_to_metadata_and_write_lineageless_fasta:
    input:
        metadata = rules.uk_filter_omitted_sequences.output.metadata,
        fasta = rules.uk_filter_omitted_sequences.output.fasta,
        previous_pangolin_lineages = config["uk_previous_pangolin_lineages"],
    output:
        metadata = config["output_path"] + "/1a/uk_latest.with_old_pangolin_lineages.csv",
        lineageless_fasta = config["output_path"] + "/1a/sequences_to_pangolin.fasta",
    log:
        config["output_path"] + "/logs/1a_uk_add_previous_pangolin_lineages_to_metadata_and_write_lineageless_fasta.log"
    # resources: mem_per_cpu=20000
    resources: mem_per_cpu=1000
    run:
        from Bio import SeqIO
        import csv

        logfile = open(str(log), "w")

        alignment = SeqIO.index(str(input.fasta), "fasta")

        lineage_dict = {}

        with open(str(input.previous_pangolin_lineages), 'r', newline = '') as lineages_in:
            reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:
                if row["taxon"] in lineage_dict:
                    logfile.write("%s occurs more than once in lineages input file\n" % row["taxon"])
                    continue

                lineage_dict[row["taxon"]] = {"lineage": row["lineage"], "pangoLEARN_version": row["pangoLEARN_version"], "probability": row["probability"]}


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.lineageless_fasta), 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["lineage", "pangoLEARN_version", "probability"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                lineage = ""
                pangoLEARN_version = ""
                probability = ""

                fasta_header = row["fasta_header"]

                if fasta_header in lineage_dict:
                    lineage = lineage_dict[fasta_header]["lineage"]
                    pangoLEARN_version = lineage_dict[fasta_header]["pangoLEARN_version"]
                    probability = lineage_dict[fasta_header]["probability"]
                else:
                    seqrec = alignment[fasta_header]
                    fasta_out.write(">" + seqrec.id + "\n")
                    fasta_out.write(str(seqrec.seq) + "\n")

                row["lineage"] = lineage
                row["pangoLEARN_version"] = pangoLEARN_version
                row["probability"] = probability

                writer.writerow(row)

        logfile.close()


rule uk_pangolin:
    input:
        fasta = rules.uk_add_previous_pangolin_lineages_to_metadata_and_write_lineageless_fasta.output.lineageless_fasta,
    params:
        outdir = config["output_path"] + "/1a/pangolin",
        tmpdir = config["output_path"] + "/1a/pangolin/tmp"
    output:
        lineages = config["output_path"] + "/1a/pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/1a_uk_add_sample_date.log"
    shell:
        """
        pangolin {input.fasta} \
            --outdir {params.outdir} \
            --tempdir {params.tmpdir} >> {log} 2>&1
        """


rule uk_add_new_pangolin_lineages_to_metadata:
    input:
        metadata = rules.uk_add_previous_pangolin_lineages_to_metadata_and_write_lineageless_fasta.output.metadata,
        lineages = rules.uk_pangolin.output.lineages,
    output:
        metadata = config["output_path"] + "/1a/uk_latest.add_new_pangolin_lineages.csv",
    log:
        config["output_path"] + "/logs/1a_uk_add_new_pangolin_lineages_to_metadata.log"
    run:
        import csv

        lineage_dict = {}

        logfile = open(str(log), "w")

        with open(str(input.lineages), 'r', newline = '') as lineages_in:
            reader = csv.DictReader(lineages_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:
                if row["taxon"] in lineage_dict:
                    logfile.write("%s occurs more than once in lineages input file\n" % row["taxon"])
                    continue

                lineage_dict[row["taxon"]] = {"lineage": row["lineage"], "pangoLEARN_version": row["pangoLEARN_version"], "probability": row["probability"]}


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                fasta_header = row["fasta_header"]

                if fasta_header in lineage_dict:
                    row["lineage"] = lineage_dict[fasta_header]["lineage"]
                    row["pangoLEARN_version"] = lineage_dict[fasta_header]["pangoLEARN_version"]
                    row["probability"] = lineage_dict[fasta_header]["probability"]

                writer.writerow(row)

        logfile.close()


rule uk_add_pillar_2:
    input:
        metadata = rules.uk_add_new_pangolin_lineages_to_metadata.output.metadata,
    output:
        metadata = config["output_path"] + "/1a/uk_latest.add_pillar_2.csv",
    log:
        config["output_path"] + "/logs/1a_uk_add_pillar_2.log"
    run:
        import csv

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["pillar_2"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                if row['collection_pillar'] == 2 or row['central_sample_id'][0:4] in ["ALDP", "CAMC", "MILK", "QEUH"]:
                    row["pillar_2"] = True
                else:
                    row["pillar_2"] = False

                writer.writerow(row)


rule uk_make_sequence_name:
    input:
        metadata = rules.uk_add_pillar_2.output.metadata,
    output:
        metadata = config["output_path"] + "/1a/uk_latest.add_sequence_name.csv",
    log:
        config["output_path"] + "/logs/1a_uk_make_sequence_name.log"
    run:
        import csv

        adm1a_to_country = {"UK-SCT": "Scotland",
                           "UK-WLS": "Wales",
                           "UK-ENG": "England",
                           "UK-NIR": "Northern_Ireland"}

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["sequence_name"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                country = adm1a_to_country[row['adm1']]
                id = row['central_sample_id']
                year = str(row['sample_date']).split("-")[0]
                name = country + "/" + id + "/" + year

                row["sequence_name"] = name

                writer.writerow(row)


rule uk_add_gisaid_accession:
    input:
        metadata = rules.uk_make_sequence_name.output.metadata,
        accessions_table = config["latest_uk_accessions"],
    output:
        metadata = config["output_path"] + "/1a/uk_latest.add_gisaid_accessions.csv"
    log:
        config["output_path"] + "/logs/1a_uk_add_gisaid_accession.log"
    run:
        import csv

        logfile = open(str(log), "w")

        accessions_dict = {}

        with open(str(input.accessions_table), 'r', newline = '') as acc_in:

            reader = csv.DictReader(acc_in, delimiter="\t", quotechar='\"', dialect = "unix")

            for row in reader:
                central_sample_id = row["central_sample_id"]
                run_name = row["run_name"]
                gisaid_accession = row["gisaid.accession"]

                if central_sample_id in accessions_dict:
                    if run_name in accessions_dict[central_sample_id]:
                        logfile.write(f'duplicate central_sample_id * run_name in accessions list: {central_sample_id} {run_name}\n')
                        continue
                    accessions_dict[central_sample_id][run_name] = gisaid_accession
                else:
                    accessions_dict[central_sample_id] = {run_name: gisaid_accession}


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["covv_accession_id"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                acc = ""
                if row["central_sample_id"] in accessions_dict:
                    if row["run_name"] in accessions_dict[row["central_sample_id"]]:
                        acc = accessions_dict[row["central_sample_id"]][row["run_name"]]

                row["covv_accession_id"] = acc

                writer.writerow(row)


rule uk_get_unmapped_genome_completeness:
    input:
        fasta = rules.uk_filter_on_sample_date.output.fasta,
        metadata = rules.uk_add_gisaid_accession.output.metadata,
    output:
        metadata = config["output_path"] + "/1a/uk_latest.annotate_genome_completeness.csv"
    log:
        config["output_path"] + "/logs/1a_uk_annotate_to_remove_duplicates.log"
    run:
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index(str(input.fasta), "fasta")

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["unmapped_genome_completeness"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                id = row["fasta_header"]

                if id in alignment:
                    seq = str(alignment[id].seq)
                    completeness = float(len(seq.replace("N", "")) / len(seq))
                    row["unmapped_genome_completeness"] = completeness
                    writer.writerow(row)
                else:
                    continue


rule uk_remove_duplicates_COGID_by_proportion_N:
    input:
        fasta = rules.uk_filter_on_sample_date.output.fasta,
        metadata = rules.uk_get_unmapped_genome_completeness.output.metadata
    output:
        fasta = config["output_path"] + "/1a/uk_latest.deduplicated_cog_id.fasta",
        metadata = config["output_path"] + "/1a/uk_latest.deduplicated_cog_id.csv"
    log:
        config["output_path"] + "/logs/1a_uk_filter_duplicates_bycovid.log"
    resources: mem_per_cpu=10000
    run:
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index(str(input.fasta), "fasta")

        dup_dict = {}

        with open(str(input.metadata), 'r', newline = '') as csv_in:
            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:

                fasta_header = row["fasta_header"]
                id = row["central_sample_id"]
                completeness = float(row["unmapped_genome_completeness"])

                if id in dup_dict:
                    if completeness > dup_dict[id]["completeness"]:
                        dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}
                    else:
                        continue
                else:
                    dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}


        tokeep = set()

        for k,v in dup_dict.items():
            tokeep.add(v["fasta_header"])


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.fasta), 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                fasta_header = row["fasta_header"]

                if fasta_header in tokeep:
                    writer.writerow(row)

                    seqrec = alignment[fasta_header]
                    fasta_out.write(">" + seqrec.id + "\n")
                    fasta_out.write(str(seqrec.seq) + "\n")

                else:
                    continue


rule uk_add_epi_week_and_day:
    input:
        metadata = rules.uk_remove_duplicates_COGID_by_proportion_N.output.metadata
    output:
        metadata = config["output_path"] + "/1a/uk_latest.epi_week.csv",
    log:
        config["output_path"] + "/logs/1a_uk_add_epi_week_and_day.log"
    run:
        import helper
        import csv

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + ["edin_epi_week", "edin_epi_day"], delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                date_str = row["sample_date"]
                epi_week = helper.date_string_to_epi_week(date_str)
                epi_day = helper.date_string_to_epi_day(date_str)

                row["edin_epi_week"] = epi_week
                row["edin_epi_day"] = epi_day

                writer.writerow(row)


rule uk_remove_duplicates_biosamplesourceid_by_date:
    input:
        fasta = rules.uk_remove_duplicates_COGID_by_proportion_N.output.fasta,
        metadata = rules.uk_add_epi_week_and_day.output.metadata
    output:
        fasta = config["output_path"] + "/1a/uk_latest.deduplicated_biosample_source_id.fasta",
        metadata = config["output_path"] + "/1a/uk_latest.deduplicated_biosample_source_id.csv"
    log:
        config["output_path"] + "/logs/1a_uk_filter_duplicates_by_biosample.log"
    # resources: mem_per_cpu=10000
    resources: mem_per_cpu=1000
    run:
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index(str(input.fasta), "fasta")

        dup_dict = {}
        tokeep = set()

        with open(str(input.metadata), 'r', newline = '') as csv_in:
            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:

                fasta_header = row["fasta_header"]
                id = row["biosample_source_id"]
                epi_day = int(row["edin_epi_day"])

                if id in ["None", "", None]:
                    tokeep.add(fasta_header)
                    continue

                if id in dup_dict:
                    if epi_day < dup_dict[id]["epi_day"]:
                        dup_dict[id] = {"fasta_header": fasta_header, "epi_day": epi_day}
                    else:
                        continue
                else:
                    dup_dict[id] = {"fasta_header": fasta_header, "epi_day": epi_day}


        for k,v in dup_dict.items():
            tokeep.add(v["fasta_header"])


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.fasta), 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                fasta_header = row["fasta_header"]

                if fasta_header in tokeep:
                    writer.writerow(row)

                    seqrec = alignment[fasta_header]
                    fasta_out.write(">" + seqrec.id + "\n")
                    fasta_out.write(str(seqrec.seq) + "\n")

                else:
                    continue


rule uk_remove_duplicates_root_biosample_by_gaps:
    input:
        fasta = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        metadata = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.metadata
    output:
        fasta = config["output_path"] + "/1a/uk_latest.deduplicated_rootbiosample.fasta",
        metadata = config["output_path"] + "/1a/uk_latest.deduplicated_rootbiosample.csv",
    log:
        config["output_path"] + "/logs/1a_uk_filter_duplicates_root_biosample_by_gaps.log"
    resources: mem_per_cpu=20000
    run:
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index(str(input.fasta), "fasta")

        dup_dict = {}
        tokeep = set()

        with open(str(input.metadata), 'r', newline = '') as csv_in:
            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:

                fasta_header = row["fasta_header"]
                id = row["root_biosample_source_id"]
                completeness = float(row["unmapped_genome_completeness"])

                if id in ["None", "", None]:
                    tokeep.add(fasta_header)
                    continue

                if id in dup_dict:
                    if completeness > dup_dict[id]["completeness"]:
                        dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}
                    else:
                        continue
                else:
                    dup_dict[id] = {"fasta_header": fasta_header, "completeness": completeness}


        for k,v in dup_dict.items():
            tokeep.add(v["fasta_header"])


        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.fasta), 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                fasta_header = row["fasta_header"]

                if fasta_header in tokeep:
                    writer.writerow(row)

                    seqrec = alignment[fasta_header]
                    fasta_out.write(">" + seqrec.id + "\n")
                    fasta_out.write(str(seqrec.seq) + "\n")

                else:
                    continue


rule uk_change_United_Kingdom_to_UK:
    input:
        metadata = rules.uk_remove_duplicates_root_biosample_by_gaps.output.metadata
    output:
        metadata = config["output_path"] + "/1a/uk_latest.united_kingdom_to_uk.csv"
    log:
        config["output_path"] + "/logs/1a_uk_change_United_Kingdom_to_UK.log"
    run:
        import csv

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                row["adm0"] = row["adm0"].replace("United Kingdom", "UK")
                writer.writerow(row)


rule uk_minimap2_to_reference:
    input:
        fasta = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/1a/alignment.sam"
    log:
        config["output_path"] + "/logs/1a_uk_minimap2_to_reference.log"
    threads: 1
    resources: mem_per_cpu=1000
    shell:
        """
        minimap2 -t {threads} -a -x asm5 {input.reference} {input.fasta} > {output.sam} 2> {log}
        """


rule uk_get_variants:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
        genbank_anno = config["reference_genbank_annotation"],
    output:
        variants = config["output_path"] + "/1a/uk.variants.csv",
    log:
        config["output_path"] + "/logs/1a_uk_get_variants.log"
    threads: 1
    shell:
        """
        /cephfs/covid/bham/raccoon-dog/programs/gofasta/gofasta sam variants -t {threads} \
            --samfile {input.sam} \
            --reference {input.reference} \
            --genbank {input.genbank_anno} \
            --outfile {output.variants} &>> {log}
        """


rule uk_get_indels:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
    output:
        insertions = config["output_path"] + "/1a/uk.insertions.txt",
        deletions = config["output_path"] + "/1a/uk.deletions.txt",
    log:
        config["output_path"] + "/logs/1a_uk_get_indels.log"
    threads: 1
    shell:
        """
        /cephfs/covid/bham/climb-covid19-jacksonb/programs/indels/indels {input.sam} {output.insertions} {output.deletions} 2> {log}
        """


rule uk_alignment:
    input:
        sam = rules.uk_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    output:
        fasta = config["output_path"] + "/1a/alignment.fasta",
    log:
        config["output_path"] + "/logs/1a_uk_alignment.log"
    threads: 1
    shell:
        """
        /cephfs/covid/bham/raccoon-dog/programs/gofasta/gofasta sam toMultiAlign -t {threads} \
            --samfile {input.sam} \
            --reference {input.reference} \
            -o {output.fasta} &>> {log}
        """


rule uk_get_snps:
    input:
        alignment = rules.uk_alignment.output.fasta,
        reference = config["reference_fasta"],
    output:
        snps = config["output_path"] + "/1a/uk.snps.csv",
    log:
        config["output_path"] + "/logs/1a_uk_get_snps.log"
    threads: 1
    shell:
        """
        /cephfs/covid/bham/climb-covid19-jacksonb/programs/snps/snps {input.reference} {input.alignment} > {output.snps} 2> {log}
        """


rule uk_type_AAs:
    input:
        metadata = rules.uk_change_United_Kingdom_to_UK.output.metadata,
        fasta = rules.uk_alignment.output.fasta,
        AAs = config["AAs"]
    output:
        metadata = config["output_path"] + "/1a/AA_finder.csv",
    log:
        config["output_path"] + "/logs/1a_uk_type_AAs.log"
    run:
        from Bio import SeqIO
        import csv

        def parse_AA_file(file):
            """
            input is in the format:
            start (1-based)
            e.g.:
            D614G,1605

            ls is a list of length-2 tuples with the format (name, position)
            position is the 1-based starting position of the codon in Wuhan-Hu-1 coordinates

            it has the same number of entries as lines in file
            """

            ls = []

            with open(file, 'r') as f:
                for line in f:
                    l = line.rstrip().split(",")
                    name, pos = l

                    ls = ls + [(name, int(pos))]

            return(ls)

        alignment = SeqIO.index(str(input.fasta), "fasta")

        AAs = parse_AA_file(str(input.AAs))

        new_columns = [x[0] for x in AAs]

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + new_columns, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                id = row["fasta_header"]
                seq = alignment[id].seq

                for entry in AAs:
                    pos = entry[1]

                    try:
                        QUERY_allele = seq[pos - 1: pos + 2].translate()
                    except:
                        QUERY_allele = 'X'

                    row[entry[0]] = QUERY_allele

                writer.writerow(row)


rule uk_type_dels:
    input:
        metadata = rules.uk_type_AAs.output.metadata,
        fasta = rules.uk_alignment.output.fasta,
        dels = config["dels"]
    output:
        metadata = config["output_path"] + "/1a/del_finder.csv",
    log:
        config["output_path"] + "/logs/1a_uk_type_dels.log"
    run:
        from Bio import SeqIO
        import csv

        def parse_del_file(file):
            """
            input is in the format:
            start (1-based), length of deletion
            e.g.:
            1605,3

            l is a list of length-3 tuples with the format (position, length, ref_allele)

            it has the same number of entries as lines in file
            """

            ls = []

            WuhanHu1 = SeqIO.read(workflow.basedir + '/MN908947.fa', 'fasta')

            with open(file, 'r') as f:
                for line in f:
                    l = line.rstrip().split(',')
                    pos, length = l
                    ref_allele = str(WuhanHu1.seq).upper()[int(pos) - 1: int(pos) - 1 + int(length)]

                    ls = ls + [(int(pos), int(length), ref_allele)]

            return(ls)

        alignment = SeqIO.index(str(input.fasta), "fasta")

        dels = parse_del_file(str(input.dels))

        new_columns = ["del_" + str(x[0]) + "_" + str(x[1]) for x in dels]

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + new_columns, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                id = row["fasta_header"]
                seq = alignment[id].seq

                for entry in dels:
                    pos = entry[0]
                    length = entry[1]
                    ref_allele = entry[2]

                    column_name = "del_" + str(pos) + "_" + str(length)

                    if seq[pos - 1: pos - 1 + length] == '-' * length:
                        genotype = 'del'
                    elif seq[pos - 1: pos - 1 + length] == ref_allele:
                        genotype = 'ref'
                    else:
                        genotype = 'X'

                    row[column_name] = genotype

                writer.writerow(row)


rule uk_filter_low_coverage_sequences:
    input:
        fasta = rules.uk_alignment.output.fasta,
        metadata = rules.uk_type_dels.output.metadata,
    params:
        min_covg = config["min_covg"],
    output:
        fasta = config["output_path"] + "/1a/alignment.low_covg_filtered.fasta",
        metadata = config["output_path"] + "/1a/alignment.low_covg_filtered.csv",
    log:
        config["output_path"] + "/logs/1a_uk_filter_low_coverage_sequences.log"
    run:
        from Bio import SeqIO
        import csv

        alignment = SeqIO.index(str(input.fasta), "fasta")

        with open(str(input.metadata), 'r', newline = '') as csv_in, \
            open(str(output.metadata), 'w', newline = '') as csv_out, \
            open(str(output.fasta), 'w') as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
            writer.writeheader()

            for row in reader:
                id = row["fasta_header"]

                if id in alignment:
                    seq = str(alignment[id].seq)
                    mapped_completeness = float(len(seq.replace("N", "")) / len(seq))
                    if mapped_completeness >= float(params.min_covg / 100):
                        writer.writerow(row)
                        fasta_out.write(">" + id + "\n")
                        fasta_out.write(seq + "\n")
                else:
                    continue


rule uk_mask_1:
    input:
        fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        mask = config["uk_mask_file"]
    output:
        fasta = config["output_path"] + "/1a/alignment.masked.fasta",
    log:
        config["output_path"] + "/logs/1a_uk_mask_1.log"
    run:
        from Bio import SeqIO

        def parse_mask_file(file):
            """
            input is in the format:
            start (1-based), mask character, regex-format string to match record.id
            e.g.:
            13402,?,^Belgium/
            d is a dictionary with the regex strings as keys and position,
            mask character and compiled regular expression as values.
            it has the same number of entries as lines in file
            """

            d = {}

            with open(file, 'r') as f:
                for line in f:
                    l = line.rstrip().split(',')
                    pos, mask_char, regex = l

                    d[regex] = {'pos': int(pos),
                                'mask_char': mask_char,
                                'regex': re.compile(regex)}

            return(d)

        mask_info = parse_mask_file(str(input.mask))

        with open(str(input.fasta), "r") as fasta_in, \
            open(str(output.fasta), "w") as fasta_out:

            for record in SeqIO.parse(fasta_in, 'fasta'):
                ID = record.id
                seq = str(record.seq)

                for entry in mask_info:
                    regex = mask_info[entry]['regex']

                    if re.search(regex, ID):
                        pos = mask_info[entry]['pos']
                        mask_char = mask_info[entry]['mask_char']

                        seq = seq[:pos - 1] + mask_char + seq[pos:]

                fasta_out.write('>' + ID + '\n')
                fasta_out.write(seq + '\n')


rule uk_trim_alignment:
    input:
        fasta = rules.uk_mask_1.output.fasta,
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
    output:
        fasta = config["output_path"] + "/1a/alignment.trimmed.fasta"
    log:
        config["output_path"] + "/logs/1a_uk_trim_alignment.log"
    run:
        from Bio import SeqIO

        strt = int(params.trim_start)
        stp = int(params.trim_end)

        with open(str(input.fasta), "r") as fasta_in, \
            open(str(output.fasta), "w") as fasta_out:

            for record in SeqIO.parse(fasta_in, "fasta"):
                seq = str(record.seq).upper()
                new_seq = ("N" * strt) + seq[strt:stp] + ("N" * (len(seq) - stp))

                fasta_out.write(">" + record.id + "\n")
                fasta_out.write(new_seq + "\n")


rule summarize_preprocess_uk:
    input:
        raw_fasta = config["latest_majora_matched_fasta"],
        fasta_filtered_by_sample_date = rules.uk_filter_on_sample_date.output.fasta,
        removed_omitted_fasta = rules.uk_filter_omitted_sequences.output.fasta,
        deduplicated_fasta_by_cogid = rules.uk_remove_duplicates_COGID_by_proportion_N.output.fasta,
        deduplicated_fasta_by_biosampleid = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        deduplicated_fasta_by_rootbiosample = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
        removed_low_covg_fasta = rules.uk_filter_low_coverage_sequences.output.fasta,
        full_unmasked_alignment = rules.uk_alignment.output.fasta,
        trimed_masked_alignment = rules.uk_trim_alignment.output.fasta,
        full_metadata = rules.uk_type_dels.output.metadata,
        final_metadata = rules.uk_filter_low_coverage_sequences.output.metadata,
        variants = rules.uk_get_variants.output.variants,
        indels = rules.uk_get_indels.output.insertions,
        snps = rules.uk_get_snps.output.snps,
    params:
        json_path = config["json_path"],
        coverage = config["min_covg"],
        time_window = config["time_window"],
    log:
        config["output_path"] + "/logs/1a_summarize_preprocess_uk.log"
    shell:
        """
        echo "> Number of sequences in raw UK fasta: $(cat {input.raw_fasta} | grep ">" | wc -l)\\n" &> {log}
        echo "> Number of sequences after removing those sampled more than {params.time_window} days ago: $(cat {input.fasta_filtered_by_sample_date} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after removing those in omissions file: $(cat {input.removed_omitted_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication by central_sample_id: $(cat {input.deduplicated_fasta_by_cogid} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication by biosample_source_id: $(cat {input.deduplicated_fasta_by_biosampleid} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after deduplication by root_source_biosample_id: $(cat {input.deduplicated_fasta_by_rootbiosample} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences after mapping and filtering on {params.coverage}% coverage: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo ">\\n" >> {log}

        mkdir -p {params.json_path}
        echo '{{"text":"' > {params.json_path}/1a_data.json
        echo "*DATAPIPE: COG-UK processing complete*\\n" >> {params.json_path}/1a_data.json
        cat {log} >> {params.json_path}/1a_data.json
        echo '"}}' >> {params.json_path}/1a_data.json
        """
        # curl -X POST -H "Content-type: application/json" -d @{params.json_path}/1a_data.json {params.grapevine_webhook}

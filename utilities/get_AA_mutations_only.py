import sys

infile = sys.argv[1]
outfile = sys.argv[2]

first = True
with open(infile, "r") as f_in, open(outfile, "w") as f_out:
    for line in f_in:
        if first:
            f_out.write(line)
            first = False
            continue

        sequence_name, sample_date, variants = line.rstrip().split(",")

        country = sequence_name.split("/")[0]
        if country not in ["England", "Wales", "Scotland", "Northern_Ireland"]:
            continue

        f_out.write(sequence_name + "," + sample_date + ",")

        variants_temp = []

        for v in variants.split("|"):
            if v.startswith("synSNP"):
                continue

            variants_temp.append(v)

        f_out.write("|".join(variants_temp) + "\n")

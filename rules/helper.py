from itertools import chain
from epiweeks import Week,Year
from datetime import datetime
from Bio import SeqIO

def date_string_to_epi_week(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi week which is cumulative total epidemiological
    weeks since 2019-12-22. Week beginning 2019-12-22 is week 0
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week:
    week = Week.fromdate(date)
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    elif week.year == 2019:
        return("0")
    else:
        cum_epi_week = week.week + len(list(chain(*[[x for x in Year(y).iterweeks()] for y in range(2020, week.year)])))
        return str(cum_epi_week)

def date_string_to_epi_day(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi day which is cumulative total days since 2019-12-22
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week week:
    week = Week.fromdate(date)
    # this is day 1 of epi-week 0:
    day_one = datetime.strptime("2019-12-22", '%Y-%m-%d').date()
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    else:
        cum_epi_day = (date - day_one).days + 1
        return str(cum_epi_day)

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


def parse_del_file(file, hu1_file):
    """
    input is in the format:
    start (1-based), length of deletion
    e.g.:
    1605,3

    l is a list of length-3 tuples with the format (position, length, ref_allele)

    it has the same number of entries as lines in file
    """

    ls = []

    WuhanHu1 = SeqIO.read(os.path.dirname(os.path.realpath(__file__)) + '/../MN908947.fa', 'fasta')

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(',')
            pos, length = l
            ref_allele = str(WuhanHu1.seq).upper()[int(pos) - 1: int(pos) - 1 + int(length)]

            ls = ls + [(int(pos), int(length), ref_allele)]

    return(ls)

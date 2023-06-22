def clean_lng(lng):
    if 'NA' in lng:
        for i, clade in enumerate(lng):
            if clade == 'NA':
                lng[i] = lng[i - 1]
    
    return lng

olines = ''
for line in open('../tsvs/tax_etoldb_up.tsv', 'r'):
    if 'code\t' not in line and line != '':
        line = line.strip().split('\t')
        line = clean_lng(line)
        olines += ('%s\td__%s;p__%s;c__%s;o__%s;f__%s;g__%s;s__%s\n' %
                   (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7]))

ofile = open('../lngs/lng_etoldb_up.txt', 'w')
ofile.write(olines)
ofile.close()

olines = ''
for line in open('../tsvs/tax_gtdb.tsv', 'r'):
    if 'code\t' not in line and line != '':
        line = line.strip().split('\t')
        line = clean_lng(line)
        olines += ('%s\td__%s;p__%s;c__%s;o__%s;f__%s;g__%s;s__%s\n' %
                   (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7]))

ofile = open('../lngs/lng_gtdb.txt', 'w')
ofile.write(olines)
ofile.close()

olines = ''
for line in open('../tsvs/tax_rvdb.tsv', 'r'):
    if 'code\t' not in line and line != '':
        line = line.strip().split('\t')

        if line[1] != 'Viruses':
            line[1] = 'Viruses'

        line = clean_lng(line)
        olines += ('%s\td__%s;p__%s;c__%s;o__%s;f__%s;g__%s;s__%s\n' %
                   (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7]))

ofile = open('../lngs/lng_rvdb.txt', 'w')
ofile.write(olines)
ofile.close()
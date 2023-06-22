lineagedict = {}
for l in open('lng_all.txt', 'r'):
    l = l.replace('\'', '')
    sid = l.split('\t')[0]
    lin = l.split('\t')[-1].replace('\n', '')
    lin = lin.replace('k__', 'd__')

    if '__NA;' in lin or lin[-1] == '__NA':
        ranks = lin.split(';')
        for i, rank in enumerate(ranks):
            if rank[-2:] == 'NA':
                ranks[i] = rank.replace('NA', '') + ranks[i - 1].rsplit('_', 1)[1]
        lineagedict[sid] = ';'.join(ranks)
    else:
        lineagedict[sid] = lin

for item in lineagedict:
    print(item + '\t' + lineagedict[item])

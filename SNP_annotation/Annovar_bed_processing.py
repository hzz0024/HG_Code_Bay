"""
infile = haplotig bed file 
outfile = bed file for regions without haplotig masking

Note input file should not have header 
Important: the last coordinate in each chromosome need further edits, 
           for example, add 
           "NC_035789.1 32278115    32650044"
           after line 
           "NC_035789.1 32128115    32278115    0   haplotig"
           The length of each chromosome is estimted by bash scrip "Annovar_cnt_length.sh", the estimated length minus 1 is the coordinate for bed file
"""
infile = 'haplotigs.bed'
outfile = 'non_haplotigs.bed'
HEADER = False
delimiter = '\t'

pre_chr = 'NULL'
pre_pos = 0
with open(infile, 'r') as f, open(outfile, 'w') as w:
    #HEADER = True
    for l in f.readlines():
        if HEADER:
            HEADER = False
            continue
        ss = l.strip().split('\t')
        chr = ss[0]
        if ss[1] == '0': # first pos is 0, skip
            pre_pos = ss[2]
            pre_chr = chr
            continue
        #ss[1] = str(int(ss[1]) - 1) # haplogit range [a, b) 
        print(chr)
        if pre_chr != chr and pre_chr != 'NULL': # not the first one, end of a chrom
            outline = pre_chr + delimiter + pre_pos + delimiter + '\n'
            w.write(outline)
            outline = chr + delimiter + '0' + delimiter + ss[1] + '\n'
            w.write(outline)
        elif pre_chr == 'NULL': # first line
            outline = chr + delimiter + '0' + delimiter + ss[1] + '\n'
            w.write(outline)
        else:
            outline = chr + delimiter + pre_pos + delimiter + ss[1] + '\n'
            w.write(outline)
        pre_chr = chr
        pre_pos = ss[2] 
        

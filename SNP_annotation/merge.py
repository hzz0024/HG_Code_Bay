infile = 'inversions.masked.bed'
outfile = 'inversions.masked.bed.merge'

hmap = {}
for i in range(1,11):
    hmap[str(i)] = []

HEADER = True
with open(infile, 'r') as f:
    for l in f:
        if HEADER:
            HEADER = False
            continue
        ss = l.strip().split('\t')
        hmap[ss[0]].append([int(ss[1]), int(ss[2])])
    
with open(outfile, 'w') as w:
    for chr in range(1,11):
        print(chr)
        arr = hmap[str(chr)]

        arr.sort(key = lambda x: x[0])
        # array to hold the merged intervals 
        m = [] 
        s = -10000
        max = -100000
        for i in range(len(arr)): 
            a = arr[i] 
            if a[0] > max: 
                if i != 0: 
                    m.append([s,max]) 
                max = a[1] 
                s = a[0] 
            else: 
                if a[1] >= max: 
                    max = a[1] 
          
        #'max' value gives the last point of  
        # that particular interval 
        # 's' gives the starting point of that interval 
        # 'm' array contains the list of all merged intervals 

        if max != -100000 and [s, max] not in m: 
            m.append([s, max]) 
        for i in range(len(m)): 
            w.write(str(chr) + '\t' + str(m[i][0]) + '\t' + str(m[i][1]) + '\n') 

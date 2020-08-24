# this python script is created to identify the common shared outliers among different SGS results.

ls1 = [] # to be compared

for i in range(1,5):
    fname = 'p_values_local_' + str(i) + '.txt'
    with open(fname, 'r') as f:
        if not ls1:
            ls1 = f.readlines()
            continue
        else:
            ls2 = f.readlines()
            print(len(ls1))
            common = []
            for l in ls1:
                for _l in ls2:
                    if _l == l:
                        common.append(_l)

            #print(len(common))
            ls1 = common

print((common))
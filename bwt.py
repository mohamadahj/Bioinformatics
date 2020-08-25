from functools import partial


# def bwt_standard(str):
#     str = str + '$'
#     table = sorted(str[i:] + str[:i] for i in range(len(str)))
#     last_column = [row[-1:] for row in table]
#     return "".join(last_column)


def radix_sort(values, key, step=0):
    if len(values) < 2:
        for value in values:
            yield value
        return

    bins = {}
    for value in values:
        bins.setdefault(key(value, step), []).append(value)

    for k in sorted(bins.keys()):
        for r in radix_sort(bins[k], key, step + 1):
            yield r


def bw_key(text, value, step):
    return text[(value + step) % len(text)]


def bwt_custom(text):
    return ''.join(text[i - 1] for i in radix_sort(range(len(text)), partial(bw_key, text)))


def C1(str, c):
    c = c.lower()
    str = "".join(sorted(str))
    i = str.find(c)
    return i


def OCC(S, c, q):
    S = S[:q]
    counter = S.count(c)
    return counter


def backward_search(S, P):
    p = len(P)
    length = len(S)
    i = p - 1
    c = P[p - 1]
    First = C1(S, c)
    Last = First + OCC(S, c, length) - 1
    # print(c, First, Last, p)
    while (First <= Last) and (i >= 1):
        c = P[i - 1]
        # print(C1(S, c))
        First = C1(S, c) + OCC(S, c, First - 1) + 1
        Last = C1(S, c) + OCC(S, c, Last + 1)
        # print(c, First, Last, OCC(S, c, First-1))
        i = i - 1
    if Last < First or First < 0:
        return 0
    else:
        return (Last - First + 1)


# str = input("Please enter the String: ")
f = open("EcoliGenome.fa", "r")
str = ""
search=""
lines2 = f.readlines()
j = 1

while(j<len(lines2)):

     #str.append((lines2[j][0:-1].lower()))
     str += lines2[j][0:-1].lower()
     j+=1

str += lines2[len(lines2)-1][-1:].lower()
str += '$'
f.close()


transformed = bwt_custom(str)

f2 = open("ReadSet5.1.fastq", "r")
lines = f2.readlines()
k = 1

with open("Result.txt","w+") as result:
	result.write("Source String:"+str+'\n')
	result.write("BWT:" + transformed+'\n')

count = 0
while(k<len(lines)):
    search = lines[k][0:-1]
    search = search.lower()
    #print(search)
    if(backward_search(transformed, search)>0):
        #print("Found.")
        with open("Result.txt", "a+") as result:
            result.write(lines[k - 1])
            result.write(search + '\n')
            result.write("Found."+'\n')
        count += 1
    #else:
        #print("Not Found.")
        # with open("Result.txt", "a+") as result:
        #     result.write("Not Found."+'\n')
    k += 4
f2.close()
# with open("Result.txt", "a+") as result:
#     result.write("Rate of Founded Strings:"+count/k)
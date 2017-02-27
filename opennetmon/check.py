rec_file = open("2015_05_09_22:20:24.recive", "r")
path_file = open("id_path.txt", "r")
path = {}
for l in path_file:
    l = l.split('->')
    if(len(l) < 2):
        continue
    path[int(l[0])] = map(int, l[1].split())
path_file.close()

tongji = {}
for k in path:
    tongji[k] = 0;
tongji[9999] = 0;
for l in rec_file:
    l = eval(l)
    tongji[int(l[0][1])] += 1
for k in tongji:
    if(tongji[k] == 0):
        print k, path[k]
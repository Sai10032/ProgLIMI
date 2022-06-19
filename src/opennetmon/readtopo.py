#!/bin/bash/env python
#-*- coding: utf-8 -*-

def readtopo():
    """
    解析result.txt文件，读入topo信息
    返回：monitors, SDN_node, paths, links
    """
    import os
    from collections import defaultdict
    from collections import namedtuple

    links = []
    monitors =[]
    SDN_node = []
    paths = []
    f = open(os.path.dirname(__file__)+'/'+'result.txt', 'r')
    line = f.readline()
    while line:
        # print(line),

        if line.find("monitors:") != -1:
            line = f.readline()
            x = line.split()
            x = [int(p)+1 for p in x]
            for i in x:
                monitors.append(i)
            # print(monitors)
            line = f.readline()
        elif line.find("SDN node:") != -1:
            line = f.readline()
            x = line.split()
            x = [int(p)+1 for p in x]
            for i in x:
                SDN_node.append(i)
            # print(SDN_node)
            line = f.readline()

        elif line.find("paths:") != -1:
            line = f.readline()
            while line:
                # print(line),
                x = line.split()
                x = [int(p)+1 for p in x]
                paths.append(tuple(x))
                line = f.readline()
            line = f.readline()

        elif line.find("links:") != -1:
            line = f.readline()
            while line.find("paths:") == -1:
                x = line.split()
                x = [int(p)+1 for p in x]
                if (x[0], x[1]) not in links and (x[1], x[0]) not in links:
                    links.append((x[0], x[1]))
                line = f.readline()
    f.close() 
    # print(paths.keys())
    return (monitors, SDN_node, paths, links)

def dijkstra(monitors, SDN_node, paths, links):
    # 使用dijkstra在SDN_node之间构建最短路径树
    pre_path = {}
    from collections import defaultdict
    dis = defaultdict(lambda:float("inf"))
    s = SDN_node[0]
    S = [s]
    dis[s] = 0
    pre_path[s] = (s, s)
    def update(i):
        for e in paths:
            if e[0] not in SDN_node or e[-1] not in SDN_node:
                #如果e的一边不是SDN_node，忽略e
                continue
            if e[0] == i:
                head = e[-1]
            elif e[-1] == i:
                head = e[0]
            else:
                continue
            if len(e)-1 + dis[i] < dis[head]:
                dis[head] = len(e)-1 + dis[i]
                if e[0] == i:
                    pre_path[head] = e
                else:
                    e = e[::-1]
                    pre_path[head] = e
    def findmin():
        min_d = float("inf")
        for k in SDN_node:
            if k in S:
                continue
            if dis[k] < min_d:
                min_d = dis[k]
                a = k
        return a
    update(s)
    while len(S) != len(SDN_node):
        i = findmin()
        S.append(i)
        update(i)

    pre = {}
    for k in pre_path:
        pre[k] = pre_path[k][0]

    # adj_path包括该节点要转发数据包的path
    adj_path = defaultdict(lambda:[])
    for p in paths:
        if p[0] in SDN_node and p[-1] not in SDN_node:
            # SDN节点到非SDN节点
            adj_path[p[0]].append(p)
        elif p[-1] in SDN_node and p[0] not in SDN_node:
            # 非SDN节点到SDN节点,反向
            p = p[::-1]
            adj_path[p[0]].append(p)
        elif p[0] in SDN_node and p[-1] in SDN_node:
            # SDN节点到SDN节点
            # print p
            if pre[p[0]] == p[-1]  or pre[p[-1]] == p[0]:
                # 存在近亲关系
                # 儿子加入到父亲的adj_path
                if p == pre_path[p[-1]]:
                    adj_path[p[0]].append(p)
                else:
                    x = p[::-1]
                    if x == pre_path[x[-1]]:
                        adj_path[x[0]].append(x)
                    else:
                        adj_path[p[0]].append(p)
            else:
                # 不存在近亲关系
                adj_path[p[0]].append(p)

    return (pre_path, dis, adj_path, pre)

if __name__ == '__main__':
    monitors, SDN_node, paths, links = readtopo()
    pre_path, dis, adj_path, pre = dijkstra(monitors, SDN_node, paths, links)
    # print adj_path
    # for k in adj_path.keys():
    #     print adj_path[k]
    print pre_path
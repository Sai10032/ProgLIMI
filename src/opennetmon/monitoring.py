#!/bin/bash/env python
#-*- coding: utf-8 -*-

"""
Listen NewMonitor evnent，install monitor flow table
"""

from pox.core import core
import pox.openflow.libopenflow_01 as of
from pox.lib.revent import *
import pox.lib.util as util
from pox.lib.recoco import Timer
from datetime import datetime
from collections import defaultdict
from collections import namedtuple
import pox.lib.packet as pkt
#from pox.openflow.of_json import flow_stats_to_list
import struct
from pox.lib.addresses import IPAddr,EthAddr
import time

log = core.getLogger()
switches = {}
switch_ports = {}
adj = defaultdict(lambda:defaultdict(lambda:None))

def _install_monitoring_path(monitor, monitors, SDN_node, paths, links, pre_path, dis, adj_path, pre):
    switches_ip, path_id, switches_mac = {}, {}, {}

    class Path_id():
        """
        随机为每一个path分配一个固定udp端口号
        """
        def __init__(self):
            self.id_dict = {}
            # 去除9999
            self.id_list = (9999,)

        def __getitem__(self, path):
            if path in self.id_dict.keys():
                return self.id_dict[path]
            x = path
            if not isinstance(path, tuple):
                x = tuple(path)
            i = hash(x)%(1<<15)
            while i in self.id_list:
                i = (i + 1)%(1<<15)
            self.id_list = self.id_list + (i,)
            self.id_dict[path] = i
            return i

        def keys(self):
            return self.id_dict.keys()
            
    path_id = Path_id()
    def topo_conf():
        """
        Each switch is configed with a ip address. You can config it with files yourself .
        """
        for k in switches.keys():
            switches_ip[k] = IPAddr((192<<24)+int(k))
            switches_mac[k] = EthAddr("aa"+ "%010d"%(k))

    def install_SDN_path():
        """
        配置全SDN topo的流表
        """
        for node in SDN_node:
            msg = of.ofp_flow_mod(match = monitor.match.clone())
            msg.match.nw_dst = switches_ip[node]
            if pre[node] == node: 
                msg.match.in_port = adj[node][monitor.dpid]
            else:
                msg.match.in_port = adj[node][pre_path[node][-2]]
                temp_path = pre_path[node]
                if temp_path[-1] in SDN_node:
                    x = temp_path[-1]
                    while pre_path[x][0] != x:
                        #log.debug(temp_path)
                        a = pre_path[x][::-1]
                        temp_path = temp_path + a[1:]
                        x = temp_path[-1]
                x = temp_path[0]
                while pre_path[x][0] != x:
                    #log.debug(temp_path)
                    temp_path = pre_path[x][:-1] + temp_path
                    x = temp_path[0]
                msg.match.tp_src = path_id[temp_path]

            #send back
            # msg.actions.append(of.ofp_action_tp_port.set_src(tp_port = path_id[pre_path[node]]))
            msg.actions.append(of.ofp_action_nw_addr.set_dst(nw_addr = switches_ip[pre[node]]))
            msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[pre_path[node][-2]]))
            msg.actions.append(of.ofp_action_output(port = of.OFPP_IN_PORT))

            # print node, msg.match.in_port
            for path in adj_path[node]:
                if path[-1] in SDN_node:
                    if len(path) == 3 and path[0] == path[-1]:
                        install_NSDN2NSDN_path(path[1:-1], path[-1], path[-1], flag = 1)
                    else:
                        install_NSDN2NSDN_path(path[1:-1], path[-1], path[-1])
                else:
                    install_NSDN2NSDN_path(path[1:], path[-1], monitor.dpid)

                temp_path = path
                if temp_path[-1] in SDN_node:
                    x = temp_path[-1]
                    while pre_path[x][0] != x:
                        #log.debug(temp_path)
                        a = pre_path[x][::-1]
                        temp_path = temp_path + a[1:]
                        x = temp_path[-1]
                x = temp_path[0]
                while pre_path[x][0] != x:
                    #log.debug(temp_path)
                    temp_path = pre_path[x][:-1] + temp_path
                    x = temp_path[0]
                # log.debug("*"*20)
                msg.actions.append(of.ofp_action_tp_port.set_src(tp_port = path_id[temp_path]))
                msg.actions.append(of.ofp_action_nw_addr.set_dst(nw_addr = switches_ip[path[-1]]))
                msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[path[-2]]))
                if adj[node][path[1]] == msg.match.in_port:
                    msg.actions.append(of.ofp_action_output(port = of.OFPP_IN_PORT))
                else:
                    msg.actions.append(of.ofp_action_output(port = adj[node][path[1]]))
            switches[node].connection.send(msg)
            # time.sleep(100000)

            msg = of.ofp_flow_mod(match = monitor.match.clone())
            msg.match.tp_src, msg.match.tp_dst, msg.match.in_port =(None,)*3
            msg.priority -= 10
            msg.match.nw_dst = switches_ip[node]
            if pre[node] != node:
                msg.actions = [of.ofp_action_nw_addr.set_dst(nw_addr = switches_ip[pre[node]])]
                msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[pre_path[node][-2]]))
                msg.actions.append(of.ofp_action_output(port = adj[node][pre_path[node][-2]]))
                back_path = pre_path[node][::-1]
                install_NSDN2NSDN_path(back_path[1:-1], back_path[-1], back_path[-1])
            else:
                msg.actions = [of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[monitor.dpid])]
                msg.actions.append(of.ofp_action_output(port = adj[node][monitor.dpid]))
            switches[node].connection.send(msg)

            msg = of.ofp_flow_mod(match = monitor.match.clone())
            msg.match.tp_src, msg.match.tp_dst, msg.match.in_port =(None,)*3
            msg.priority -= 1
            msg.match.nw_dst = switches_ip[node]
            if pre[node] != node:
                msg.match.in_port = adj[node][pre_path[node][-2]]
                msg.actions = [of.ofp_action_nw_addr.set_dst(nw_addr = switches_ip[pre[node]])]
                msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[pre_path[node][-2]]))
                msg.actions.append(of.ofp_action_output(port = of.OFPP_IN_PORT))
                switches[node].connection.send(msg)

    def install_NSDN2NSDN_path(path, dst, next_hop, flag = None):
        """
        Install flow table without SDN node in "path" to send packet to "dst",
        when the packet arrive at path[-1], next_hop decide packet to which port in path[-1]
        """
        if len(path) == 0:
            return
        msg = of.ofp_flow_mod(match = monitor.match.clone())
        msg.match.nw_dst = switches_ip[dst]
        msg.match.tp_src, msg.match.tp_dst, msg.match.in_port = None, None, None
        # print msg.match
        # time.sleep(1000)
        for i in range(len(path)):
            msg.actions = []
            if i+1 <= len(path)-1:
                msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[path[i+1]]))
                msg.actions.append(of.ofp_action_output(port = adj[path[i]][path[i+1]]))
                switches[path[i]].connection.send(msg)
            else:
                msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[next_hop]))
                if flag == 1:
                    msg.priority += 10
                    msg.match.in_port = adj[path[i]][next_hop]
                    msg.actions.append(of.ofp_action_output(port = of.OFPP_IN_PORT))
                else:
                    msg.actions.append(of.ofp_action_output(port = adj[path[i]][next_hop]))
                switches[path[i]].connection.send(msg)
            # if(path[i] == 12):
            #     print msg.actions[-1].port, "flag:%s"%flag
    def install_r_flow():
        """
        Using a r SDN node  send and receive probe packet
        """
        msg = of.ofp_flow_mod()
        msg.match = monitor.match.clone()
        for node in SDN_node:
            if pre[node] == node:
                break
        msg.actions.append(of.ofp_action_nw_addr.set_dst(nw_addr = switches_ip[node]))
        msg.actions.append(of.ofp_action_output(port = adj[monitor.dpid][node]))
        for k in monitors:
            if k in SDN_node:
                continue
            for p in paths:
                if p[0] == k:
                    msg.actions.append(of.ofp_action_nw_addr.set_dst(nw_addr = switches_ip[p[-1]]))
                    msg.actions.append(of.ofp_action_tp_port.set_src(tp_port = path_id[p]))
                    msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = switches_mac[p[0]]))
                    msg.actions.append(of.ofp_action_output(port = adj[monitor.dpid][k]))
        switches[monitor.dpid].connection.send(msg)

        msg = of.ofp_flow_mod()
        msg.match = monitor.match.clone()
        msg.match.in_port, msg.match.nw_dst, msg.match.tp_dst, msg.match.tp_src = (None,)*4
        msg.match.nw_src = monitor.ip
        msg.actions.append(of.ofp_action_nw_addr.set_dst(nw_addr = monitor.ip))
        msg.actions.append(of.ofp_action_nw_addr.set_src(nw_addr = IPAddr('10.0.0.2')))
        msg.actions.append(of.ofp_action_dl_addr.set_dst(dl_addr = monitor.hw))
        msg.actions.append(of.ofp_action_dl_addr.set_src(dl_addr = EthAddr("00:00:00:00:00:02")))
        msg.actions.append(of.ofp_action_output(port = monitor.switch_port))
        switches[monitor.dpid].connection.send(msg)

    topo_conf()

    #配置连接所有发送和接受probe数据包的ovs，这样看上在测量topo看上去就是多个host在发送和接收probe packet
    install_r_flow()

    for path in paths:
        if path[0] not in SDN_node and path[-1] not in SDN_node: #paths without SND node
            install_NSDN2NSDN_path(path, path[-1], monitor.dpid)
    install_SDN_path()

    with open("id_path.txt", "w") as f:
        # 将每条路径的udp端口号和路径写入文件，udphandler.py根据这个区分不同测量路径
        f.write("id_path:\n")
        for p in path_id.keys():
            f.write("%6s -> "%path_id[p])
            f.write("%s "*len(p)%p)
            f.write("\n")
    with open("topo.txt", "w") as f:
        for k in adj:
            for j in adj[k]:
                f.write("%s %s %s\n"%(k, j, adj[k][j]))
def _build_monitoring_topo(monitor):

    from readtopo import readtopo, dijkstra
    monitors, SDN_node, paths, links = readtopo()
    (pre_path, dis, adj_path, pre) = dijkstra(monitors, SDN_node, paths, links)
    _install_monitoring_path(monitor, monitors, SDN_node, paths, links, pre_path, dis, adj_path, pre)

class Monitoring (object):        
    def __init__ (self):
        log.debug("Monitoring coming up")
        def startup():
            # 监听_handle_LinkEvent事件
            core.openflow_discovery.addListeners(self)
            # 监听_handle_ConnectionUp事件
            core.openflow.addListeners(self)
            # 监听_handle_NewMonitor事件
            core.opennetmon_handle_PacketIn.addListeners(self)
            log.debug("Monitoring started")
        core.call_when_ready(startup, 'opennetmon_handle_PacketIn') #Wait for opennetmon-forwarding to be started

    def _handle_ConnectionUp (self, event):
        switches[event.connection.dpid] = event
      
    def _handle_NewMonitor(self, event):
        log.debug("_handle_NewMonitor")
        _build_monitoring_topo(event)
    
    def _handle_LinkEvent(self, event):
        link = event.link
        if event.added:
            adj[link.dpid1][link.dpid2] = link.port1
            
def launch ():
    # 注册Monitoring组件
    core.registerNew(Monitoring)
    

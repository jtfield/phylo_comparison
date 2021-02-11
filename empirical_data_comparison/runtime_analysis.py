#! /usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict
import numpy as np
import pandas as pd
import dendropy
from dendropy.calculate import treecompare
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gon_phy_times')
    parser.add_argument('--ep_times')
    parser.add_argument('--snippy_times')
    return parser.parse_args()

def ep_gon_process_times(time):
    total_time_per_run_in_seconds = 0
    time = time[0]
    split_time = time.split('m')
    minutes = int(split_time[0])
    minutes_in_seconds = minutes * 60
    seconds = int(split_time[1])

    total_time_per_run_in_seconds = minutes_in_seconds + seconds

    return total_time_per_run_in_seconds

def process_times(time_file):
    snippy_id = "Walltime"
    snippy_id_compile = re.compile(snippy_id)
    ep_and_gon_id = "real"
    ep_gon_id_compile = re.compile(ep_and_gon_id)
    
    ep_gon_time = "\d+m\d+"
    ep_gon_time_compile = re.compile(ep_gon_time)
    snippy_minutes = "\d+\sminutes"
    snippy_seconds = "\d+\sseconds"
    snip_min_compile = re.compile(snippy_minutes)
    snip_sec_compile = re.compile(snippy_seconds)

    ep_gon_time_in_seconds_list = []
    snippy_time_in_seconds_list = []

    line_count = 0
    with open(time_file, 'r') as times:
        for line in times:
            line_count+=1
            ep_gon_find = re.findall(ep_gon_id_compile, line)
            snippy_find = re.findall(snippy_id_compile, line)
            if ep_gon_find:
                # print(ep_gon_find)
                ep_gon_time_find = re.findall(ep_gon_time_compile, line)
                if ep_gon_time_find:
                    # print(ep_gon_time_find)
                    time_per_run = ep_gon_process_times(ep_gon_time_find)
                    ep_gon_time_in_seconds_list.append(time_per_run)
            
            
            elif snippy_find:
                run_time = []
                run_time_in_seconds = 0
                # print(snippy_find)
                snip_sec_find = re.findall(snip_sec_compile, line)
                if snip_sec_find:
                    # print(snip_sec_find)
                    run_time.append(int(snip_sec_find[0].strip(' seconds')))
                    snip_min_find = re.findall(snip_min_compile, line)
                    if snip_min_find:
                        # print(snip_min_find)
                        run_time.append(int(snip_min_find[0].strip(' minutes')) * 60)
                if len(run_time) > 1:
                    run_time_in_seconds = run_time[0] + run_time[1]
                    snippy_time_in_seconds_list.append(run_time_in_seconds)
                elif len(run_time) == 1:
                    run_time_in_seconds = run_time[0]
                    snippy_time_in_seconds_list.append(run_time_in_seconds)
                
    if len(ep_gon_time_in_seconds_list) > 0:
        return ep_gon_time_in_seconds_list
    elif len(snippy_time_in_seconds_list) > 0:
        return snippy_time_in_seconds_list
    
    
def fig_gen(ep_times, gon_phy_times, snippy_times):

    del ep_times[-1]
    del ep_times[-1]
    del ep_times[-1]

    ep_columns = ["time"]
    snip_columns = ["time"]
    gon_phy_columns = ["time"]
    ep_df = pd.DataFrame(ep_times, columns=ep_columns)
    gon_phy_df = pd.DataFrame(gon_phy_times, columns=gon_phy_columns)
    snippy_df = pd.DataFrame(snippy_times, columns=snip_columns)

    ep_df['method'] = 'Extensiphy'
    gon_phy_df['method'] = 'De novo'
    snippy_df['method'] = 'Snippy'

    frames = [ep_df, snippy_df, gon_phy_df]

    combined_df = pd.concat(frames)
    # print(combined_df)
    combined_df.reset_index(drop=True, inplace=True)
    # print(combined_df)
    combined_df['Minutes'] = combined_df['time'].div(60)
    # print(combined_df)

    # ax = sns.violinplot(x='method', y="minutes", data=combined_df, inner=None, color=".8")
    ax = sns.stripplot(x="method", y="Minutes", data=combined_df, linewidth=1, jitter=0.3)
    ax.set_title( "Time for Individual Sequence Assembly" , size = 18 )
    ax.set_xlabel( "Method" , size = 12 )
    ax.set_ylabel( "Minutes per Assembly" , size = 12 ) 
    plt.show()

    # ep_avg = ep_df['ep_time'].mean()
    # snip_avg = snippy_df['snip_time'].mean()
    # gon_avg = gon_phy_df['gon_phy_time'].mean()

    # print(ep_avg / 60)
    # print(gon_avg / 60)
    # print(snip_avg / 60)

    # ep_df['time'].div(60)
    # snip_min_df = snippy_df['time'].div(60)
    # gon_min_df = gon_phy_df['time'].div(60)

    # print(ep_df)

    # combined_df = pd.DataFrame([ep_min_df['ep_time'], snip_min_df['snip_time'], gon_min_df['gon_phy_time']]).transpose()
    # combined_df.dropna(axis=0, how='any', inplace=True)
    
    # ax = sns.stripplot(x=combined_df.index, y=combined_df.values)
    
    # plt.show()


    # fig=plt.figure()
    # ax=fig.add_axes([0,0,1,1])
    # ax.scatter(combined_df['gon_phy_time'].max(), combined_df['ep_time'], color='r')
    # ax.scatter(combined_df['gon_phy_time'].max(), combined_df['snip_time'], color='b')
    # ax.scatter(combined_df['gon_phy_time'].max(), combined_df['gon_phy_time'], color='y')
    # ax.set_xlabel('time to assemble')
    # ax.set_ylabel('number of assemblies')
    # ax.set_title('scatter plot')

    # plt.show()

    # print(combined_df)

    # n_bins = 60

    # fig, axs = plt.subplots(1, 3, tight_layout=True, sharey=True, figsize=(10,10))
    # fig.suptitle('Individual Sequences Time-to-assemble', fontsize=18, fontweight='bold')
    # plt.ylim(0,400)

    # axs[0].hist(ep_min_df['times_in_seconds'], bins=n_bins, label="A")
    # axs[0].set_title('Extensiphy', fontsize=16)

    # axs[1].hist(gon_min_df['times_in_seconds'], bins=n_bins, label="B")
    # axs[1].set_title('De novo assembly', fontsize=16)

    # axs[2].hist(snip_min_df['times_in_seconds'], bins=n_bins, label="C")
    # axs[2].set_title('Snippy', fontsize=16)

    # fig.text(0.5, 0.01, 'Minutes Per Assembly', ha='center', va='center', fontsize=14)
    # fig.text(0.01, 0.5, 'Number of Sequences', ha='center', va='center', rotation='vertical', fontsize=14)

    # plt.show()

def main():
    args = parse_args()

    ep_times = process_times(args.ep_times)
    # print(ep_times)

    snippy_times = process_times(args.snippy_times)
    # print(snippy_times)

    gon_phy_times = process_times(args.gon_phy_times)
    # print(gon_phy_times)

    fig_gen(ep_times, gon_phy_times, snippy_times)





if __name__ == '__main__':
    main()
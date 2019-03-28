import os
import csv
import pandas as pd


#retry_li = ['1moq','4cox', '1lpz', '1vbm','1rv1','1jak']


with open('info.csv') as csv_file:
    pdb_info = csv.reader(csv_file, delimiter=',')
    all_li = []
    for line in pdb_info:
        pdb = line[0].lower()
        # if pdb in retry_li:
        #     continue
        lig_id = line[3].split('_')[-1].lower()
        if line[2] == 'v' and os.path.exists('pdb_files/{}/{}/tractability_single/scores.csv'.format(pdb[1:3], pdb)):
            target_df = pd.read_csv('pdb_files/{}/{}/tractability_single/scores.csv'.format(pdb[1:3], pdb))
            print(target_df)
            target_df['pdb'] = pdb
            target_df['tractability'] = line[1]
            all_li.append(target_df)

    all_df = pd.concat(all_li)

    all_df.to_csv('scores_single.csv')

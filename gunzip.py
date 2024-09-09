import os
import gzip

# 指定待解压的目录
target_dir = '/dssg/home/acct-clslhz/clslhz/hws/Alpha-GEMs/struct_data/yeast'
# 指定解压后的目录
output_dir = '/dssg/home/acct-clslhz/clslhz/hws/Alpha-GEMs/struct_data/yeast'

# 遍历目录下的所有gz文件，并解压到指定目录
for file_name in os.listdir(target_dir):
    if file_name.endswith('.gz') or file_name.endswith('.cif'):
        os.remove(target_dir+'/'+file_name)
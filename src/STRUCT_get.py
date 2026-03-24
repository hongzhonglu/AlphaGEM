# import os
# import pandas as pd
# from selenium import webdriver
# from selenium.webdriver.common.by import By
# from selenium.webdriver.chrome.service import Service
# from selenium.webdriver.chrome.options import Options
# import time
#
# # 初始化全局变量
# name = 'rhto'
# input_path = f'working/{name}/{name}.xlsx'
# output_dir = f'./struct_data/taryeast/{name}'
# os.makedirs(output_dir, exist_ok=True)  # 确保输出目录存在
#
# # 设置 Chrome 浏览器驱动路径
# chrome_driver_path = '/usr/bin/chromedriver'  # 修改为实际的 chromedriver 路径
#
# # 设置 Chrome 浏览器选项
# chrome_options = Options()
# chrome_options.add_argument("--headless")  # 无头模式，后台运行，不弹出浏览器窗口
# chrome_options.add_argument("--disable-gpu")  # 禁用 GPU
# chrome_options.add_argument("--no-sandbox")  # 无沙盒模式
#
# # 启动 Chrome 浏览器
# service = Service(chrome_driver_path)
# driver = webdriver.Chrome(service=service, options=chrome_options)
#
# # 打开 AlphaFold 网站
# driver.get("https://www.alphafold.ebi.ac.uk/")
#
# # 等待页面加载
# time.sleep(3)
#
# # 读取 Excel 文件
# taryeast = pd.read_excel(input_path)
#
#
# # 下载函数
# def download(gene_id, url, output_path):
#     try:
#         # 打开目标页面
#         driver.get(url)
#         time.sleep(0.25)  # 等待页面加载
#
#         # 等待文件下载完成
#         time.sleep(0.25)  # 根据文件大小调整等待时间
#
#         print(f"Downloaded: {gene_id}")
#         return 1
#     except Exception as e:
#         print(f"Error downloading {gene_id}: {e}")
#         return 0
#
#
# # 多线程下载（模拟多线程下载操作）
# def main():
#     futures = []
#
#     for i, row in taryeast.iterrows():
#         gene_id = row[0]
#         if isinstance(row[2], float):  # 跳过无效数据
#             continue
#
#         url = f"https://alphafold.ebi.ac.uk/files/AF-{gene_id}-F1-model_v4.pdb"
#         output_path = f"{output_dir}/AF-{gene_id}-F1-model_v4.pdb"
#
#         # 如果文件已存在，跳过下载
#         if os.path.exists(output_path):
#             print(f"File already exists: {output_path}")
#             continue
#
#         # 下载文件
#         download(gene_id, url, output_path)
#
#     # 清理无效文件
#     for file_name in os.listdir(output_dir):
#         if file_name.endswith('.1'):
#             os.remove(os.path.join(output_dir, file_name))
#             print(f"Deleted: {file_name}")
#
#
# # 执行下载任务
# if __name__ == "__main__":
#     main()
#
# # 关闭浏览器
# driver.quit()

import os
import pandas as pd
import requests  # 需要安装: pip install requests
import time

# 初始化全局变量
name = 'pichi'
input_path = f'working/{name}/{name}.xlsx'
output_dir = f'./struct_data/taryeast/{name}'
os.makedirs(output_dir, exist_ok=True)

# 读取 Excel
taryeast = pd.read_excel(input_path)


def download_file(url, save_path):
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # 检查请求是否成功

        # 写入文件
        with open(save_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Downloaded: {save_path}")
        return True
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False


def main():
    for i, row in taryeast.iterrows():
        gene_id = row[0]
        # 根据实际情况调整 row 的索引
        if pd.isna(gene_id):
            continue

        url = f"https://alphafold.ebi.ac.uk/files/AF-{gene_id}-F1-model_v6.pdb"
        output_path = os.path.join(output_dir, f"AF-{gene_id}-F1-model_v4.pdb")

        # 检查文件是否存在
        if os.path.exists(output_path):
            print(f"File already exists: {output_path}")
            continue

        # 下载
        download_file(url, output_path)

        # 礼貌性延时，防止被封IP
        time.sleep(0.5)


if __name__ == "__main__":
    main()

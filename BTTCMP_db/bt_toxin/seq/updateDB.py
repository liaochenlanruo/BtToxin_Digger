def get_unique_headers(file_path):
    """读取文件中以'>'开头的行，返回'>'后面内容的集合"""
    headers = set()
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 提取'>'后面的内容（包括可能的空格和其他字符）
                header = line[1:]
                headers.add(header)
    return headers

# 输入文件路径
file1 = 'bt_toxin20251104.fas'
file2 = 'all_app_cry_cyt_gpp_mcf_mpf_mpp_mtx_pra_prb_spp_tpp_txp_vip_vpa_vpb_xpp_fasta_sequences.txt'
output_file = 'unique_headers.txt'

# 获取两个文件中的header集合
headers1 = get_unique_headers(file1)
headers2 = get_unique_headers(file2)

# 计算各自独有的header
unique_to_file1 = headers1 - headers2
unique_to_file2 = headers2 - headers1

# 写入输出文件
with open(output_file, 'w') as out_f:
    out_f.write(f"### Unique headers in {file1} ###\n")
    for header in sorted(unique_to_file1):
        out_f.write(f">{header}\n")
    
    out_f.write(f"\n### Unique headers in {file2} ###\n")
    for header in sorted(unique_to_file2):
        out_f.write(f">{header}\n")

print(f"处理完成，结果已保存至 {output_file}")
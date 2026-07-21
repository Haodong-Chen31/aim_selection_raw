#!/usr/bin/env python3

import pandas as pd
import argparse
from collections import defaultdict
import sys

def filter_snps_by_fst_and_ld(fst_file, ld_file, fst_threshold=0.01, ld_threshold=0.1, output_file="filtered_snps.txt"):
    """
    根据FST值和LD关系筛选SNP
    
    参数:
    fst_file: FST文件路径，格式为 染色体 位置 FST值
    ld_file: LD文件路径，PLINK格式，包含SNP对和r²值
    fst_threshold: FST阈值，默认0.01
    ld_threshold: LD阈值，默认0.1
    output_file: 输出文件路径
    """
    
    # 1. 读取FST数据
    print("读取FST数据...")
    try:
        fst_df = pd.read_csv(fst_file, sep='\s+', header=0, 
                            names=['CHR', 'POS', 'WEIR_AND_COCKERHAM_FST'])
    except Exception as e:
        print(f"错误: 无法读取FST文件 {fst_file}: {e}")
        sys.exit(1)
    
    # 创建SNP标识符
    fst_df['SNP_ID'] = fst_df['CHR'].astype(str) + ':' + fst_df['POS'].astype(str)
    
    # 筛选FST > threshold的SNP
    high_fst_snps = fst_df[fst_df['WEIR_AND_COCKERHAM_FST'] > fst_threshold].copy()
    print(f"FST > {fst_threshold} 的SNP数量: {len(high_fst_snps)}")
    
    if len(high_fst_snps) == 0:
        print("警告: 没有找到FST值超过阈值的SNP")
        # 创建空的结果文件
        pd.DataFrame(columns=['CHR', 'POS', 'WEIR_AND_COCKERHAM_FST']).to_csv(output_file, sep='\t', index=False, header=False)
        return
    
    # 2. 读取LD数据
    print("读取LD数据...")
    try:
        ld_df = pd.read_csv(ld_file, sep='\s+', header=0,
			   names=['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2'])
    except Exception as e:
        print(f"错误: 无法读取LD文件 {ld_file}: {e}")
        sys.exit(1)
    
    # 标准化列名（适应不同PLINK版本）
    ld_columns = {}
    for col in ld_df.columns:
        if 'CHR' in col and '_A' in col:
            ld_columns[col] = 'CHR_A'
        elif 'CHR' in col and '_B' in col:
            ld_columns[col] = 'CHR_B'
        elif 'SNP' in col and '_A' in col:
            ld_columns[col] = 'SNP_A'
        elif 'SNP' in col and '_B' in col:
            ld_columns[col] = 'SNP_B'
        elif 'BP' in col and '_A' in col:
            ld_columns[col] = 'BP_A'
        elif 'BP' in col and '_B' in col:
            ld_columns[col] = 'BP_B'
        elif 'R2' in col:
            ld_columns[col] = 'R2'
    
    ld_df = ld_df.rename(columns=ld_columns)
    
    # 筛选LD > threshold的SNP对
    high_ld_pairs = ld_df[ld_df['R2'] > ld_threshold].copy()
    print(f"LD > {ld_threshold} 的SNP对数量: {len(high_ld_pairs)}")
    
    # 3. 构建SNP网络图
    print("构建SNP关联网络...")
    snp_to_fst = dict(zip(high_fst_snps['SNP_ID'], high_fst_snps['WEIR_AND_COCKERHAM_FST']))
    
    # 创建邻接表
    snp_graph = defaultdict(list)
    for _, row in high_ld_pairs.iterrows():
        snp_a = f"{row['CHR_A']}:{row['BP_A']}"
        snp_b = f"{row['CHR_B']}:{row['BP_B']}"
        
        if snp_a in snp_to_fst and snp_b in snp_to_fst:
            snp_graph[snp_a].append(snp_b)
            snp_graph[snp_b].append(snp_a)
    
    # 4. 贪心算法选择SNP
    print("执行贪心算法选择SNP...")
    selected_snps = set()
    visited = set()
    
    # 按FST值降序排序
    sorted_snps = sorted(snp_to_fst.keys(), 
                        key=lambda x: snp_to_fst[x], reverse=True)
    
    for snp in sorted_snps:
        if snp in visited:
            continue
            
        # 选择当前SNP（FST最高的未访问SNP）
        selected_snps.add(snp)
        visited.add(snp)
        
        # 标记所有高LD的邻居为已访问
        for neighbor in snp_graph[snp]:
            visited.add(neighbor)
    
    print(f"最终选择的SNP数量: {len(selected_snps)}")
    
    # 5. 输出结果
    result_df = high_fst_snps[high_fst_snps['SNP_ID'].isin(selected_snps)]
    result_df = result_df[['CHR', 'POS', 'WEIR_AND_COCKERHAM_FST']].sort_values(['CHR', 'POS'])
    
    # 保存结果
    result_df.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"结果已保存到: {output_file}")
    
    # 显示统计信息
    print(f"\n筛选统计:")
    print(f"总SNP数量: {len(result_df)}")
    print(f"染色体分布:")
    print(result_df['CHR'].value_counts().sort_index())
    
    return result_df

def main():
    """主函数：解析命令行参数并执行筛选"""
    parser = argparse.ArgumentParser(description='根据FST值和LD关系筛选SNP')
    parser.add_argument('--fst', '-f', required=True, help='FST文件路径（格式：染色体 位置 FST值）')
    parser.add_argument('--ld', '-l', required=True, help='LD文件路径（PLINK格式）')
    parser.add_argument('--fst-threshold', '-t', type=float, default=0.01, 
                       help='FST阈值（默认：0.01）')
    parser.add_argument('--ld-threshold', '-r', type=float, default=0.1,
                       help='LD阈值（默认：0.1）')
    parser.add_argument('--output', '-o', default='filtered_snps.txt',
                       help='输出文件路径（默认：filtered_snps.txt）')
    
    args = parser.parse_args()
    
    print("=" * 50)
    print("SNP筛选工具")
    print("=" * 50)
    print(f"FST文件: {args.fst}")
    print(f"LD文件: {args.ld}")
    print(f"FST阈值: {args.fst_threshold}")
    print(f"LD阈值: {args.ld_threshold}")
    print(f"输出文件: {args.output}")
    print("=" * 50)
    
    # 执行筛选
    filter_snps_by_fst_and_ld(
        fst_file=args.fst,
        ld_file=args.ld,
        fst_threshold=args.fst_threshold,
        ld_threshold=args.ld_threshold,
        output_file=args.output
    )

if __name__ == "__main__":
    main()

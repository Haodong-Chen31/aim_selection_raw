# coding:utf-8
# @Time : 2022/03/08 17:45
# @File : pcaRploter
# @Version: v4.5
# @Author : zky

import argparse
import sys
import os


def joinn(var, sep=None, chr=""):
    if sep:
        return chr.join(var.split(sep=sep))
    else:
        return chr.join(var.split())


Usage = """
this procedure will create a Rscript file for further pca plotting
use -h for detailed
poplist format:
    ====pop1====
    ind1
    ind2
    ====pop2====
    ind3
    ind4
"""
print(Usage)
parser = argparse.ArgumentParser()
group = parser.add_argument_group()
group.add_argument('-p', '--poplist', help='specify poplist, default=poplist.txt',
                   default='poplist.txt', type=str)
group.add_argument('-o', '--output', help='specify output file, default=smartpca.r',
                   default='smartpca.r', type=str)
group.add_argument('-i', '--input', help='specify input plot file, default=plot.txt',
                   default='plot.txt', type=str)
group.add_argument('-pdf', '--pdfout', help='specify PDF output file, default=smartpca.pdf',
                   default='smartpca.pdf', type=str)
group.add_argument('-legend', '--legend', help='specify legend PDF output file, default=legend.pdf',
                   default='legend.pdf', type=str)
args = parser.parse_args()
print('poplist :', args.poplist)
print('output : ', args.output)
print('input plot file : ', args.input)
print('PDF output : ', args.pdfout)
print('Legend PDF output : ', args.legend)

# 读取poplist原始内容
with open(args.poplist, 'r') as f:
    poptext_original = f.readlines()

# 读取输入文件的第一列（人群信息）
try:
    with open(args.input, 'r') as f:
        input_lines = f.readlines()
        # 提取输入文件中所有的人群标签
        input_populations = set()
        for line in input_lines:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 1:
                    input_populations.add(parts[0])
except FileNotFoundError:
    print(f"错误: 输入文件 '{args.input}' 不存在!")
    sys.exit(1)

print("\n" + "="*50)
print("清理poplist文件...")
print("="*50)

# 创建一个清理后的poplist
cleaned_poplist = []
current_group_header = None
current_group_individuals = []

# 处理poplist原始内容
for line in poptext_original:
    line = line.rstrip('\n')
    
    if line.startswith('===='):
        # 保存前一个组（如果有的话）
        if current_group_header is not None:
            if current_group_individuals:
                # 如果组中有个体，保存整个组
                cleaned_poplist.append(current_group_header)
                cleaned_poplist.extend(current_group_individuals)
            else:
                print(f"删除空组: {current_group_header.strip('=')}")
        
        # 开始新组
        current_group_header = line
        current_group_individuals = []
    elif line.startswith('#') or not line.strip():
        # 跳过注释行和空行
        continue
    else:
        # 检查个体是否在输入文件中存在
        if line in input_populations:
            current_group_individuals.append(line)
        else:
            print(f"删除不存在于输入文件中的个体: {line}")

# 处理最后一组
if current_group_header is not None:
    if current_group_individuals:
        cleaned_poplist.append(current_group_header)
        cleaned_poplist.extend(current_group_individuals)
    else:
        print(f"删除空组: {current_group_header.strip('=')}")

# 如果没有有效的数据，退出
if not cleaned_poplist:
    print("错误: 清理后poplist为空!")
    print("请检查poplist文件和输入文件的人群对应关系")
    sys.exit(1)

# 保存清理后的poplist到临时文件
import tempfile
temp_poplist = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_cleaned.txt')
for line in cleaned_poplist:
    temp_poplist.write(line + '\n')
temp_poplist.close()

print(f"清理后的poplist已保存到: {temp_poplist.name}")
print(f"有效组数: {sum(1 for line in cleaned_poplist if line.startswith('===='))}")
print(f"有效个体数: {sum(1 for line in cleaned_poplist if not line.startswith('===='))}")

# 使用清理后的poplist重新解析
poptext = cleaned_poplist
pop_list = []
j = -1
for i in poptext:
    if i.startswith('='):
        j += 1
        pop_list.append([])
        pop_list[j].append(i.strip("="))
    else:
        pop_list[j].append(i)

# 生成R脚本
pca = open(args.output, 'w')
pca.writelines(f"""pdf("{args.pdfout}")
data=read.table("{args.input}")
pop=data[,1]
ind=data[,2]
PC1=data[,3]
PC2=data[,4]
# layout(matrix(c(1,2)),widths=c(2,1),heights=c(3,1))
# par(mar=c(4,4,1,5))
plot(PC1,PC2,type="n")
frame=data.frame(PC1=PC1,PC2=PC2,name=pop,region=pop)
light=c("#FFDAB9", "#ABF4DC", "#D3D3D3", "#b1c8f3", "#ffbbea", "#c8b0f5", "#e2ffa8", "#c1fffa", "#f3f3a9", "#fa9c9c")
deep=c("#CCCC00", "#3787B4", "#00CD66", "#6B8E23", "#D02090", "#FF8C69", "#905aff", "#20B2AA", "#FFC534", "#FFFF00", "#6F96F2", "#8D5223", "#FF1493")


# COLOR SETTINGS
""")

light = 1
deep = 1
for i in range(0, len(pop_list)):
    # 所有人群组都使用非实心点，不区分target和非target
    if "Ancient" in pop_list[i][0]:
        pca.writelines(
            joinn(joinn(pop_list[i][0], "-")) + f'_col=deep[{str(deep)}]\n')
        deep += 1
    else:
        pca.writelines(
            joinn(joinn(pop_list[i][0], "-")) + f'_col=light[{str(light)}]\n')
        light += 1

# PLOT POPs
pca.writelines("\n# PLOT\n")
for i in pop_list:
    pca.writelines(f'# {i[0]}\n')
    for j in range(1, len(i)):
        pca.writelines(f'reg="{i[j]}"\n')
        if j < 26:
            col_tmp = joinn(joinn(i[0], "-"))
            # 使用pch从1开始的非实心符号，bg=NA表示不填充
            pca.writelines(
                'points(subset(frame,region==reg)$PC1,subset(frame,region==reg)$PC2,pch={},col={}_col,bg=NA,cex=0.6)\n'.format(
                    str(j), col_tmp))
        else:
            col_tmp = joinn(joinn(i[0], "-"))
            pca.writelines(
                'points(subset(frame,region==reg)$PC1,subset(frame,region==reg)$PC2,pch={},col={}_col,bg=NA,cex=0.6)\n'.format(
                    str(j + 7), col_tmp))
pca.writelines("dev.off()\n\n# LEGEND\n")
pca.writelines(f'pdf("{args.legend}")\n')
pca.writelines('par(mar=c(5,5,5,5))\nplot.new()\n')
# pops
pops = []
for i in pop_list:
    pops.append('","'.join(i))
pops = '"'+'","'.join(pops)+'"'
pca.writelines(f'pops=c({pops})\n')
# cols,borders,symb,fonts
cols = []
bgs = []
symb = []
fonts = []
for i in pop_list:
    cols.append('NA')
    bgs.append('NA')
    symb.append('NA')
    fonts.append('2')
    le = len(i) - 1
    pp = joinn(joinn(i[0], "-")) + "_col"
    cols.append(f"rep({pp}, {le})")
    bgs.append('NA')  # 所有点都不填充，使用bg=NA
    fonts.append(f"rep(1, {le})")
    [symb.append(str(j + 1)) if j < 25 else symb.append(str(j + 8)) for j in range(le)]
cols = ",".join(cols)
bgs = ",".join(bgs)
symb = ",".join(symb)
fonts = ",".join(fonts)

pca.writelines(f'cols=c({cols})\n')
pca.writelines(f'bgs=c({bgs})\n')
pca.writelines(f'symb=c({symb})\n')
pca.writelines(f'fonts=c({fonts})\n')
pca.writelines(
    'legend("top",pops,pch=symb,pt.bg=bgs,col=cols,ncol=3,cex=0.4,pt.cex=0.4,pt.lwd=0.7,bty="o",text.font=fonts,xpd=TRUE)\n')
pca.writelines('dev.off()\n')
pca.close()

print("\n" + "="*50)
print("清理结果总结:")
print("="*50)
print(f"原始poplist文件: {args.poplist}")
print(f"清理后poplist文件: {temp_poplist.name}")
print(f"生成的R脚本: {args.output}")

# 显示最终的人群组信息
print("\n最终使用的人群组:")
for i, pop_group in enumerate(pop_list):
    pop_name = pop_group[0]
    individual_count = len(pop_group) - 1
    print(f"  {i+1}. {pop_name}: {individual_count} 个人群")

print("\n" + "="*50)
print(f"请运行: Rscript {args.output}")
print("="*50)

# 可选：询问是否删除临时文件
#response = input("\n是否删除临时清理文件? (y/n, 默认y): ")
#if response.lower() != 'n':
#    os.unlink(temp_poplist.name)
#    print(f"已删除临时文件: {temp_poplist.name}")
#else:
#    print(f"临时文件保留在: {temp_poplist.name}")

# 自动删除临时文件
os.unlink(temp_poplist.name)
print(f"\n已自动删除临时清理文件")

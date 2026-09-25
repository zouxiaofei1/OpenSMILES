# IUPAC Red Book 2005（无机命名）PDF → Markdown 转换产物

本目录存放 IUPAC《Nomenclature of Inorganic Chemistry: IUPAC Recommendations 2005》
（Red Book）各章 PDF 及其 Markdown 转换稿。中译本见 `../inorganic_translated/`。

## 文件

| 文件 | 说明 |
|---|---|
| `IR-N.pdf` | 原书分章 PDF（IR-1 ~ IR-11） |
| `IR-N.md` | 转换稿（保留上下标、粗斜体、表格、结构式插图） |
| `Tables_appendix.md` | 原 PDF 中附在 IR-11 之后的 Table I–X 与主题索引（未翻译） |
| `images/` | 从 PDF 抽取的插图：结构式、流程图、旋转排版的表格页 |
| `tools/pdf2md.py` | 转换脚本（pymupdf） |

转换命令：

```
cd E:/dev/chem && .venv/Scripts/python docs/iupac/inorganic/tools/pdf2md.py \
    docs/iupac/inorganic/IR-6.pdf docs/iupac/inorganic/IR-6.md
```

## 转换要点

原 PDF 由 3B2 排版系统生成，字体子集的 `ToUnicode` 映射不可靠，转换脚本做了针对性修复：

- **错误字形表**：同一字体族内被错误映射的字符（如 `þ`→`+`、`ð`→`(`、`¼`→`=`、
  `Z`→`η`、`s`→`σ`、`e`→`Δ` 等）按族逐字形渲染核对后改正，见表 `FIX`。
- **上下标**：按行组主字号与基线偏移判定上标/下标，优先转 Unicode 上下标
  （`Fe₃O₄`、`SO₄²⁻`），无对应字符时退回 `<sub>`/`<sup>`（如 `V<sub>M</sub><sup>x</sup>`）。
- **词间空格**：PDF 有时丢失空格，按字符级位置间隔补回。
- **插图**：矢量结构式（键线多为文字字形，括线为矢量路径）按括线聚类、扩入邻近
  短标签行后整块渲染为 PNG，并**按页内纵向位置插入正文**；图形内部的标签文字不再
  进入正文。只由横线（表格框线）构成、且无结构键线的区域不当作插图，避免把表格
  文字误吞成图片；旋转 90° 排版的表格页则整页保留为图片。
- **表格**：按行切单元格、列起点链式聚类成 markdown 表格，跨页同列数表格自动合并。
- **行末连字符**：用语料高频词表判定是断词还是真连字符。

## 已知取舍

- **旋转 90° 排版的表格页**（IR-8 的 Table IR-8.1/IR-8.2、IR-6 的 Table IR-6.1、
  IR-11 附表）按整页插图保留，未做文字化；`IR-N.md` 中对应位置为图片引用。
- **带框线的结构式表格**（IR-10 的表 IR-10.1–10.4、IR-7 的部分表）同样以整页/整块
  图片保留，因为这些表的首栏是结构式，文字层无法还原。
- IR-9、IR-10 的立体结构式数量多（合计近 200 张），均以插图形式保留。
- 原书正文中的印刷错误照录，不擅改。

# MyGeno: interpret your 23andMe data
MyGeno is a Python module which reads and interprets 23andMe raw genotype data. It Support simple wellness report from 23andMe, and searching your genotype with related traits at SNPedia.

## Quick-look
```python
>>> from MyGeno import Geno
>>> g = Geno("/path/to/your/raw/data.txt")
```

- `Geno.report_wellness()` prints out wellness report provided by 23andMe API.

```python
>>> g.report_wellness()
```
```
=============================== WELLNESS REPORT =============================== 

< Alcohol Flush Reaction >
--------------------------------------------------------------------------------
The marker we tested comes in two different forms, the G variant and the A varia
nt. The A variant results in an enzyme that is less efficient at breaking down a
cetaldehyde. The A variant is also known as c.1510G>A, Glu487Lys, and Glu504Lys.
 This marker has been studied the most in people of East Asian descent.

  rsID: rs671
  Yours: G;G - homozygous
  Detail: 2 | Alcohol Flush: Normal, doesn't flush. Normal hangovers. Normal ris
k of Alcoholism. Normal risk of Esophageal Cancer. Disulfiram is effective for a
lcoholism. | 60.2% of JPT

...
```

- `Geno.type_id(rsid)` prints out specific SNP of input rsid and its related trait written in SNPedia.

```python
>>> from MyGeno import Geno
>>> g = Geno("/path/to/your/raw/data.txt")
```

- `Geno.report_wellness()` prints out wellness report provided by 23andMe API.

```python
>>> g.type_id("rs28897696")
```
```
rs28897696 | G;G | 0 | normal | 0.0% of JPT
```
<div>
<table border="1" class="dataframe">
<thead>
<tr style="text-align: right;">
<th></th>
<th>rsid</th>
<th>chr</th>
<th>pos</th>
<th>genotype</th>
</tr>
</thead>
<tbody>
<tr>
<th>492260</th>
<td>rs28897696</td>
<td>17</td>
<td>38469446</td>
<td>GG</td>
</tr>
</tbody>
</table>
</div>

- `Geno.type_trait(trait)` prints out all SNPs related to input trait, written in SNPedia.

```python
>>> g.type_trait("baldness")
```

```
rs6152 | G;G | 0.5 | able to go bald | 100.0% of JPT
rs2223841 | T;T | 1.2 | more likely to go bald before age 40 | 100.0% of JPT
rs6152 | G;G | 0.5 | able to go bald | 100.0% of JPT
rs2223841 | T;T | 1.2 | more likely to go bald before age 40 | 100.0% of JPT
rs2180439 | C;T | 2.5 | Increased risk of Male Pattern Baldness. | 44.2% of JPT
rs11683401 | C;T |  | --% of JPT
rs2073963 | G;G | 2.5 | increased risk of baldness | 18.6% of JPT
rs1511061 | T;T |  | --% of JPT
 * more details at https://www.snpedia.com/index.php/baldness
```
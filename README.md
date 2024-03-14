# Curve Alignment By Best Correlation

Calculate the best cross-correlation between experiment and theory for gene identification

![image](reference/schematics.png)

In the j-th block, I calculate the 1) maximum correlation, 2) the best scale factor and 3) best offset using an experimental curve and a theoretical unzipping curve of the j-th gene of the whole genome. A whole library of theoretical curves has been established using the [unzip_GPU](https://github.com/Taomihog/unzip_GPU)


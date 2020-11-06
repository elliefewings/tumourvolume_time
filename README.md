# tumourvolume_time
## R script to perform t.test and GLM analysis for time based tumour x treatment analysis

### Requirements
 * R - packages should be automatically installed but if there are errors please open R and install the packages manually
  
## Usage
```
$ ./tumourvol_time.R [options]

Options:
 -i INPUT, --input=INPUT
 Path to CSV containing formatted tumour volume data (see example) [required]

 -c CONTROL, --control=CONTROL
 Name of control treatment [required]
 
 -h, --help
 Show this help message and exit

```
#See example_input.csv for an example of how to layout the input csv file

## Output
Generated data will include:
1) plot and data for t.test of the effect of treatment vs control at each time point per tumour model
2) plot and data for t.test of the effect of treatment vs control at each time point for all tumour models combined
3) plot and data for a GLM of the effect of treatment vs control per tumour model
3) plot and data for a GLM of the effect of treatment vs control for all tumour models combined


For the GLM, effect sizes/ are in the output csv file under the column "glm.standard.coefficient"

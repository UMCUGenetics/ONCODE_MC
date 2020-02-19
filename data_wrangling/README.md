# Data wrangling

## Intro

We have data from the paper [Linnekamp paper](https://www.nature.com/articles/s41418-017-0011-5). It is gene expression measured using microarrays.

## Scripts

`process_datasets.R`: takes the microarray data and adds fake noise, bad samples, batch effects... to play with in the notebook. This fake raw data is stored in the rawdata folder.

`data_wrangling.Rmd`: the notebook of the masterclass.

## Masterclass main points

- Normalization (Marc)
- Scaling
- Bad samples (Marc)
- Data integration
- Feature selection
- Batch correction
- Missing values (Marc)

## TODO

- Finish the notebook.
- Consider other main points.
- Consider adding more noise to the data.
- Consider summarizing the probesets into genes to reduce the data.
- Keywords for Jeroen to introduce in his lecture.
- Perhaps some introductory slides.
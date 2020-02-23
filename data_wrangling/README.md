# Data wrangling

## Intro

We have data from the paper [Linnekamp paper](https://www.nature.com/articles/s41418-017-0011-5). It is gene expression measured using microarrays.

## Scripts

`download_dataset.R`: downloads the data.

`parse_datasets.R`: extracts the downloaded data and adds the CMS subtypes to it.

`process_datasets.R`: takes the microarray data and adds fake noise, bad samples, batch effects... to play with in the notebook. This fake raw data is stored in the rawdata folder.

`data_wrangling.Rmd`: the notebook of the masterclass.

## Masterclass main points

- Normalization (Marc)
  - Polish and comment on other methods
- Scaling
- Bad samples (Done)
- Data integration
- Feature selection
- Batch correction (Done)
- Missing values (Done)

## TODO

- Keywords for Jeroen to introduce in his lecture.
- Perhaps some introductory slides.
- Prepare the second dataset to be predicted.
- Prepare the functions for the SVM blackbox.
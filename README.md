## Oncode masterclass: Introduction to machine learning in cancer genomics - Practical topics

### 1. Data wrangling

Data preparation accounts for about 80% of the work of data scientists (source CrowdFlower data science report 2016). In this workshop, we will expose the common pitfalls one can do when working with cancer genomics datasets, from the treatment of missing data to normalisation and feature selection (for ML purpose).
### 2. Clustering:

Clustering basics: we will explore commonly used clustering algorithms such as hierarchical clustering and k-means clustering, and apply this knowledge to identify clusters in an expression dataset of breast cancer patients.
### 3. Random forest:

Random forest (RF) is an ensemble-based classification method that aggregates the predictive outcomes of multiple decision trees to overcome the bias and error of weak-learning classifiers. In this session, we will explore how to implement RF models in Python and apply it to genomic datasets.
### 4. Deep Learning 1:

Image classification with Deep Learning in Python: join the journey of building a simple CNN model with Keras – this example will also entail some more generic ML know-how about cross-validation, imbalanced class case handling and more.
### 5. Deep Learning 2:

Data augmentation with Generative Adversarial Networks (GANs): an introduction to GANs and how they can be applied to augment a training set for the purpose of improving the performance of a classifier
### 6. Transfer Learning:
Machine learning models are often used to predict labels on a dataset different from the one used for training, even for datasets that diverge only slightly from the training data, this usually leads to a drop in predictive performance. In this session, we will explore approaches to train and correct a predictive machine learning model in order to achieve robust performance on datasets from other domains, e.g. training a drug response predictor on cell lines and then employing it to predict patient response.


More information on the event, please visit website: https://www.oncode.nl/events/masterclass-introduction-to-machine-learning-in-cancer-genomics

### How to run the notebooks:

To run the notebooks you can choose to use:
- [Google Colab](https://colab.research.google.com/notebooks/intro.ipynb);
- [Binder](https://mybinder.org/);
- [miniconda](https://docs.conda.io/en/latest/miniconda.html) by installing the
 conda environment and run jupyter locally from your laptop:

```console
# download Miniconda3 installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# install Conda (respond by 'yes')
bash miniconda.sh
# update Conda
conda update -y conda
# create & activate new environment with installed dependencies
git clone https://github.com/UMCUGenetics/ONCODE_MC.git
cd ONCODE_MC
conda env create -n oncode-mc -f environment.yaml
conda activate oncode-mc
jupyter notebook
```

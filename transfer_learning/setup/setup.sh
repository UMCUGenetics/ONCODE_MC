conda create -n myenv python=3.7
conda activate oncode_mc
conda install pytorch torchvision torchaudio cudatoolkit=10.2 -c pytorch
pip install seaborn scikit-learn pandas
conda install ipykernel
python -m ipykernel install --user --name='Python (Oncode Masterclass)'
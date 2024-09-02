# MetAssimulo 2: a web app for simulating realistic 1D & 2D Metabolomic 1H NMR spectra 

## Setup Instructions
### 1. Clone the Repository
```bash
git clone [https://github.com/yourusername/yourrepository.git](https://github.com/yanyan5420/MetAssimulo_2.git)
cd MetAssimulo_2
```
### 2. Download the Input Data
```bash
curl -c ./cookie -s -L "[https://drive.google.com/uc?export=download&id=YOUR_FILE_ID](https://drive.google.com/file/d/12k-nbIcoME5GSkD8u1iYUcNzg6l-7msW/view?usp=share_link)" | grep -o 'confirm=[^&]*' | sed 's/confirm=//g' > ./confirm.txt
curl -Lb ./cookie "[https://drive.google.com/uc?export=download&id=YOUR_FILE_ID](https://drive.google.com/file/d/12k-nbIcoME5GSkD8u1iYUcNzg6l-7msW/view?usp=share_link)&confirm=$(<./confirm.txt)" -o Input.zip
unzip Input.zip -d ./
rm Input.zip
```

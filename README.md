# MitoEdit: a tool for targeted mitochondrial base editing

## About
**MitoEdit** is a novel Python workflow designed for targeted base editing of the mitochondrial DNA (mtDNA). This tool streamlines the selection of optimal targeting windows for base editing, helping researchers efficiently target mtDNA mutations that are linked to various diseases. By enabling targeted base editing of the mitochondrial genome, MitoEdit will speed up the study of mtDNA-related diseases, help with preclinical drug testing, and enable therapeutic approaches to correct pathogenic mutations.

![Rough workflow](imgs/MitoEdit_Pipeline.png)

## Overview
MitoEdit lets users input DNA sequences in text format, specify the target base position, and indicate the desired modification. The tool processes this information to identify candidate target windows, list the number and position of potential bystander edits, and find optimal flanking TALE sequences, where applicable. All the results are provided in a structured format including detailed logs to track progress.

### Pipelines
The current MitoEdit workflow established four different pipelines based on the editing patterns observed in the following base editor systems: 
1. [Mok2020_G1333](https://rdcu.be/dXUIG)
2. [Mok2020_G1397](https://rdcu.be/dXUIG)
3. [Mok2022_G1397_DddA11](https://rdcu.be/dXUIm)
4. [Cho_sTALEDs](https://www.sciencedirect.com/science/article/pii/S0092867422003890)

A detailed description for each pipeline can be found in the [README](pipelines/README.md) file in the *pipelines* folder.

**Note**: To use the evolved DddA6 variant from the [Mok 2022 paper](https://rdcu.be/dXUIm), you can use the `Mok2020_G1397` pipeline.

### TALE-NT Tool
The workflow uses the [TALE-NT tool](https://tale-nt.cac.cornell.edu/) to identify optimal flanking TALE arrays around each candidate target window.

## Quick Start
If you have all the necessary tools installed, you can use these commands to access MitoEdit:

#### For targeting the human mitochondrial DNA:
```
python mitocraft.py <position> <reference_base> <mutant_base>
```
**Note:** The `position` is based on the [human mitochondrial genome](https://www.ncbi.nlm.nih.gov/nuccore/251831106) sequence (NC_012920.1)

#### For targeting any other DNA sequence:
```
python mitocraft.py --input_file <input_DNA_file> <position> <reference_base> <mutant_base>
```

## Prerequisites
- Python 3.x
- Required Python packages:
  - `pandas` 
  - `openpyxl`
- Download and install [Conda](https://docs.anaconda.com/anaconda/install/) if not already installed. Follow the prompts to complete installation.

## Installation

#### 1. Clone the repository:
```
git clone https://github.com/Kundu-Lab/mitoedit.git
cd mitoedit
```
#### 2. Set up the Conda environment:
* Use the provided `environment.yml` file to create the Conda environment:
```
conda env create -f environment.yml
```
- Alternatively, you can create the conda environment manually:
```
conda create -n run_talen_env python=2.7.18 biopython=1.70
```
**Note:** The Conda environment should be named `run_talen_env` for MitoEdit to correctly use the TALE-NT tool. There is no need to activate the conda environment, as the pipeline will automatically use `run_talen_env` if it is installed in Conda. This environment is specifically for the TALE-NT Tool.

## Web Interface

MitoEdit includes a web interface that makes it easy to analyze DNA sequences through your browser. Here's how to set it up:

### System Requirements
- Python 3.10 or newer
- Conda

### Setup Instructions

1. **Install Dependencies**
```bash
# Create environments for web server and TALE-NT tool
conda env create -f web-environment.yml
conda env create -f talen-environment.yml
```

2. **Run the Web Server**
```bash
# Set the required password for web interface access
export MITOEDIT_PASSWORD=your_secure_password

# Start the web server
conda run -n mitoedit-web python web/main.py
```

The web interface will be available at http://localhost:8000, where you can:
- Input DNA sequence position and bases
- Upload custom DNA sequence files
- View and download analysis results in Excel format

### Troubleshooting

If you encounter issues:
- Check the terminal output for error messages
- Review the log file: `logging_main.log`
- Ensure both `mitoedit-web` and `run_talen_env` environments are properly created
- Verify all required files and directories exist

### Docker Installation

MitoEdit is also available as a Docker container, which provides an isolated environment with all dependencies pre-installed.

#### Prerequisites
- Docker installed on your system ([Install Docker](https://docs.docker.com/get-docker/))

#### Building and Running with Docker

1. **Build the Docker image:**
```bash
sudo docker build -t mitoedit .
```

2. **Run the container:**

For HTTP (port 80):
```bash
sudo docker run -d -p 80:80 -e MITOEDIT_PASSWORD=your_secure_password mitoedit
```

For HTTPS (port 443):
```bash
sudo docker run -d -p 443:443 -e PORT=443 -e MITOEDIT_PASSWORD=your_secure_password mitoedit
```

The web interface will be available at:
- HTTP: http://localhost:80 (when using port 80)
- HTTPS: https://localhost:443 (when using port 443)

#### Stopping/Terminating the Docker Container

To stop or terminate the running MitoEdit Docker container:

1. **List all running containers to find the container ID:**
```bash
sudo docker ps
```

2. **Stop the container:**
```bash
sudo docker stop <container_id>
```
   For example: `sudo docker stop f11e4cc1615e`

3. **Remove the container (optional, if you want to completely remove it):**
```bash
sudo docker rm <container_id>
```

#### Docker Troubleshooting
- Ensure Docker is running on your system
- Check if port 80 is available (if not, you can map to a different port using `-p 8000:80` for example)
- For HTTPS, check if port 443 is available (if not, you can map to a different port using `-p 8443:443` for example)
- If you encounter permission issues, ensure you have the necessary privileges to run Docker commands

### Running as a Systemd Service

You can set up MitoEdit as a systemd service to run automatically on system boot:

1. **Create the service file:**

The repository includes a `mitoedit.service` file that defines how to run the Docker container as a systemd service.

2. **Install the service:**

Run the provided installation script with sudo:
```bash
sudo ./install_mitoedit_service.sh
```

3. **Managing the service:**

Once installed, you can manage the service using standard systemd commands:
```bash
sudo systemctl start mitoedit.service
sudo systemctl stop mitoedit.service
sudo systemctl restart mitoedit.service
sudo systemctl status mitoedit.service
```

4. **Accessing the web interface:**

When the service is running, access the MitoEdit web interface at:
```
https://localhost:443
```

## What input parameters does MitoEdit require?
MitoEdit requires the following parameters:
### Required Data Parameters: 

#### 1. **Position**: The position of the base to target (1-based index).

#### 2. **Reference Base**: The original base at the specified position (A, T, C, or G).

#### 3. **Mutant Base**: The desired base conversion (A, T, C, or G).

### Optional Data Parameter:
#### Input File:
- The path to a file (.txt / .fasta) containing the DNA sequence. 
- If not provided, MitoEdit will use the [human mtDNA sequence](https://www.ncbi.nlm.nih.gov/nuccore/251831106) from NCBI by default.

## What does MitoEdit output?
MitoEdit generates the following outputs:

### Output Directories
MitoEdit organizes the output files in the following directories:
- fasta: Contains the FASTA file of the 60bp sequence adjacent to the target base.
- pipeline_windows: Stores target windows generated from each individual pipeline.
- all_windows: Contains a combined list of target windows from all pipelines.
- talen: Stores the output from the TALE-NT tool.
- matching_output: Includes TALE sequences for applicable target windows.
- ####  **final_output: Contains the final consolidated output file.**

**Note:** Check the `final_output` directory to see the final results. The other directories are stored under the `running` directory.

### Output Files

#### Excel Files:
- `{pipeline}_{position}.xlsx`: Lists the target windows generated from each pipeline.
- `all_windows_{position}.xlsx`: Contains all target windows combined from all pipelines.
- `matching_tales_{position}.xlsx`: Summary of optimal flanking TALE sequences for each applicable target window. 
- #### **`final_{position}.xlsx`: Consolidated final result, including the target windows, position and number of bystander edits, and optimal flanking TALE sequences where applicable.**

#### FASTA File:
- `adjacent_bases_{position}.fasta`: Contains the sequence adjacent to the target base, extending 30bp on each side.

#### TALE-NT File:
- `TALENT_{position}.txt`: Contains the output from TALE-NT Tool describing the optimal flanking TALE sequences possible.

## Usage
- For a full list of parameters, use the --help flag from the command line.
```
python mitocraft.py --help
python mitocraft.py -h
```
- To run the tool, use the following commands:
```
python mitocraft.py <position> <reference_base> <mutant_base>
python mitocraft.py --input_file <DNA.txt> <position> <reference_base> <mutant_base>
```
### Examples
#### To target the human mitochondrial DNA:

```
python mitocraft.py 11696 G A
```
**Expected Output:**
When you run this command, MitoEdit generates an Excel file named `final_11696.xlsx` in the `final_output` directory. This file includes two spreadsheets: **All_Windows** and **Bystander_Effect**, with the first five rows from each shown below.

**Note**: The [  ] represents the target base and {  } represent bystander edits.

**1. All_Windows Sheet**
| Pipeline|  Position |Reference Base | Mutant Base | Window Size | Window Sequence | Target Location | Number of bystanders | Position of Bystanders | Optimal Flanking TALEs | Flag CheckBystanderEffect |
|--------------|----------------------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|
|Mok2022_G1397_DddA11|		11696	|G|	A	|14bp|	GCA[G]TCATT{C}TCAT|	Position 4 from the 5' end	|1|	[11702]	|FALSE	| -
|Mok2022_G1397_DddA11	|	11696	|G	|A|	14bp|	CGCA[G]TCATT{C}TCA|	Position 5 from the 5' end|	1|	[11702]|	FALSE|	 -
|Mok2022_G1397_DddA11	|11696|G	|A|	14bp|	GCGCA[G]T{C}ATTCTC|	Position 6 from the 5' end|	1	|[11698]	|FALSE	| -
|Mok2022_G1397_DddA11	|	11696	|G|	A|	14bp	|GGCGCA[G]T{C}ATTCT|	Position 7 from the 5' end|	1|	[11698]|	FALSE	| -

**Note:** If the column `Flag CheckBystanderEffect=TRUE`, you should manually check the results for potential amino acid changes caused by neighbouring bystanders on the same codon.

**2. Bystanders_Effects Sheet**
|Bystander Position|	Reference Base|	Mutant Base|	Location On Genome|	Predicted Mutation Impact	|SNV Type|	AA Variant|	Functional Impact|	MutationAssessor Score|
|---------|---------|---------|	---------|	---------	|---------|---------|---------|---------|
|11698	|C	|T	|Complex 1|	Predicted Benign		|synonymous SNV|V313V|		||
|11702	|C	|T	|Complex 1|	Predicted Pathogenic		|nonsynonymous SNV|V315F|medium|3.44|	||
|11704	|C	|T	|Complex 1|	Predicted Benign		|synonymous SNV|V315L|		||


#### To target any other DNA sequence:
```
python mitocraft.py --input_file test.txt 33 G A
```

**Expected Output:**
When using an input file, the generated Excel file will contain only one spreadsheet, similar to the following: (example taken from the file provided in the [test file](test/input/test.txt) folder)

**1. All_Windows Sheet**
| Pipeline| Position |Reference_Base | Mutant Base | Window Size | Window Sequence | Target Location| Number of bystanders | Position of Bystanders | Optimal Flanking TALEs | Flag_CheckBystanderEffect |LeftTALE1 |	RightTALE1|LeftTALE2|RightTALE2 |
|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------|---------|-------------|
|Mok2020_G1397|	33|	G|	A|	14bp|	TG{G}[G]A{G}AACT{C}TCT|	Position 4 from the 5' end	|3|	[32, 35, 40]	|FALSE|	TRUE|				
|Mok2020_G1397|	33|	G|	A|	14bp|	CTG{G}[G]A{G}AACTCTC|	Position 5 from the 5' end|	2|	[32, 35]|	FALSE|	TRUE	|			
|Mok2020_G1397|	33|	G|	A|	14bp|	ACTG{G}[G]AGAACTCT|	Position 6 from the 5' end	|1|	[32]	|FALSE|	TRUE|				
|Mok2020_G1397|	33|	G|	A|	14bp|	TACTG{G}[G]AGAACTC|	Position 7 from the 5' end	|1|	[32]	|TRUE	|TRUE	|T TACCCCCCACTATTAACC	|TCTGTGCTAGTAACC A	|T ACCCCCCACTATTAACC	|TCTGTGCTAGTAACC A|

**Note**: When optimal flanking TALE sequences are found, the sequence is added to the `LeftTALE` and `RightTALE` columns respectively. The impact of bystander edits is not provided when using an input DNA file.

## Notes
- **Input File Formatting**: Ensure that your input file is correctly formatted, with the reference base matching the base at the specified position.
- **Supported File Formats**: MitoEdit can read `.fasta` and `.txt` formats for DNA sequences.
- **Minimum Upload Sequence Length**: The input file must contain at least 35 bases, covering the target base on either side, for accurate processing.
- **Output File Name Conflicts**: Check for existing output files with the same name before running MitoEdit to prevent overwriting and errors.
- **Logging:** MitoEdit logs its progress and any issues encountered during execution in `logging_main.log.`
- **Species Support**: While the tool is designed primarily for human mtDNA, other DNA sequences can also be uploaded and used.
- **Modifying TALE-NT Workflow**: If no matching flanking TALE sequences are identified, consider modifying the TALE-NT parameter by setting `FILTER = 2`. This will identify all TALE pairs targeting any base in the target window, not just those for the target base. For further information, refer to the [TALE-NT FAQs](https://tale-nt.cac.cornell.edu/faqs).


## How to Cite?
If you use MitoEdit in your research, please cite the following paper:
- [MitoEdit: a pipeline for optimizing mtDNA base editing and predicting bystander effects](https://doi.org/10.1101/2025.01.22.634390)


## Contact
If you have any questions, please do not hesitate to contact me at:
- Devansh Shah: Devansh.Shah@stjude.org, shah.16@iitj.ac.in

## License
This project is licensed under the MIT License. See the `LICENSE` file for details. You are free to use and modify it, but please give credit to the original authors.

## Contributing
We welcome contributions! If you want to help improve MitoEdit, please fork the repository and submit a pull request.

## TALE-NT License
Copyright (c) 2011-2015, Nick Booher njbooher@gmail.com and Erin Doyle edoyle@iastate.edu.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE

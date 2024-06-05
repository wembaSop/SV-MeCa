# Use the Miniconda base image
FROM condaforge/mambaforge:latest

# Set the working directory to /workspace
WORKDIR /workspace
RUN mkdir input

# Clone the GitHub repository
RUN apt-get update && \
    apt-get install -y less nano git && \
    git clone https://github.com/wembaSop/SV-MeCa.git

# deactivate warning outdated conda
RUN conda config --set notify_outdated_conda false

#move to the SV-MeCa repository
WORKDIR /workspace/SV-MeCa

#create the required environments
RUN for file in $(ls conda/*.yml);do mamba env create -f $file;done

# Install required packages
RUN conda config --add channels bioconda
RUN mamba install -y nextflow=23.10.1 singularity=3.8.6

# fix version bug of samtools
RUN mamba install -n breakdancer samtools=1.17
RUN mamba install -n pindel samtools=1.17

# Set the entrypoint to a shell and change to the workspace directory
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["cd /workspace && /bin/bash"]

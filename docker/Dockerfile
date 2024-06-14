# Use the Miniconda base image
FROM condaforge/mambaforge:latest

# Set the working directory to /workspace
WORKDIR /workspace
RUN mkdir input

# Clone the GitHub repository & install tardis dependencies
RUN apt-get update -y --fix-missing 
RUN echo 'tzdata tzdata/Areas select Etc' | debconf-set-selections
RUN echo 'tzdata tzdata/Zones/Etc select UTC' | debconf-set-selections
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata
RUN apt-get upgrade -y
RUN apt-get install -y less nano git && \
    apt-get install git make gcc zlib1g-dev libbz2-dev liblzma-dev gnuplot-nox -y && \
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

# install tardis native
RUN git clone --recursive https://github.com/BilkentCompGen/tardis.git /workspace/tardis
RUN git clone https://github.com/BilkentCompGen/mrfast.git /workspace/mrfast
WORKDIR /workspace/tardis
RUN make libs && make
WORKDIR /workspace/mrfast
RUN make
ENV PATH="/workspace/tardis:/workspace/mrfast:${PATH}"
RUN apt-get remove git -y
RUN apt-get autoremove -y
WORKDIR /workspace

# Set the entrypoint to a shell and change to the workspace directory
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["cd /workspace && /bin/bash"]

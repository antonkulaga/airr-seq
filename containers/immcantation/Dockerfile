FROM immcantation/lab:devel
USER root
RUN python3 -m pip install jupyterlab jupyterlab-lsp jupyterlab_code_formatter
RUN pip install papermill
#RUN dnf -y install 'dnf-command(copr)' && dnf -y copr enable iucar/cran && dnf -y install R-CoprManager
ADD install.R /tmp/
RUN R -f /tmp/install.R
USER magus
WORKDIR /home/magus/notebooks
COPY src/* /home/magus/notebooks/
EXPOSE 8888
ENTRYPOINT []
CMD jupyter lab .
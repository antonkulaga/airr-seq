FROM immcantation/lab:devel
USER root
RUN python3 -m pip install jupyterlab jupyterlab-lsp jupyterlab_code_formatter
USER magus
EXPOSE 8888
ENTRYPOINT []
CMD jupyter lab .
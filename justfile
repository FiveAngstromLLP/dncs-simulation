export DNCS_FOLDER := justfile_directory()

@run:
    echo "Working on directory {{DNCS_FOLDER}}"
    python src/main.py

@install:
    pip install -r requirements.txt
    pip install maturin
    maturin develop --release -m dncs/python/Cargo.toml

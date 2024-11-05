export DNCS_FOLDER := justfile_directory()

@run:
    echo "Working on directory {{DNCS_FOLDER}}"
    python src/main.py

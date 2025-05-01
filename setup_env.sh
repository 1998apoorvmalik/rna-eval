#!/bin/bash

ENV_DIR="venv"

# Check if the script is being sourced (optional warning)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "[WARNING] Don't run this script directly if you want to activate the environment. Use: source $0"
fi

# Create virtual environment if it doesn't exist
if [ ! -d "$ENV_DIR" ]; then
    echo "[INFO] Creating virtual environment in $ENV_DIR..."
    python3 -m venv "$ENV_DIR"
    
    echo "[INFO] Activating virtual environment and installing dependencies..."
    source "$ENV_DIR/bin/activate"

    pip install --upgrade pip
    pip install hopcroftkarp

    # Make main scripts executable
    chmod +x eval_sd.py
    chmod +x eval_ed.py
    chmod +x eval_seq_identity.py
    chmod +x ct2db.py
    chmod +x ct2seq.py

    echo "[INFO] Setup complete."
    echo "To activate the environment later, run: source $ENV_DIR/bin/activate"
else
    echo "[INFO] Virtual environment already exists."
fi

# Activate the virtual environment
source "$ENV_DIR/bin/activate"
echo "[INFO] Virtual environment activated."

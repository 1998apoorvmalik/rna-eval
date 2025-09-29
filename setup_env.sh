#!/bin/bash

ENV_DIR="venv"

# Create virtual environment if it doesn't exist
if [ ! -d "$ENV_DIR" ]; then
    echo "[INFO] Creating virtual environment in $ENV_DIR..."
    python3 -m venv "$ENV_DIR"
    
    echo "[INFO] Installing dependencies into virtual environment..."
    "$ENV_DIR/bin/pip" install --upgrade pip
    "$ENV_DIR/bin/pip" install hopcroftkarp

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
    echo "To activate the environment later, run: source $ENV_DIR/bin/activate"
fi

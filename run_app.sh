#!/bin/bash

# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate

# Upgrade pip and install dependencies
pip install --upgrade pip
pip install -r requirements.txt

# Run the Streamlit app
streamlit run main.py

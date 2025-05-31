from fastapi import FastAPI, UploadFile, File, Form
from fastapi.responses import JSONResponse
from SequenceAnalyzer import SequenceAnalyzer
import tempfile
import shutil

app = FastAPI()

@app.post("/analyze")
async def analyze(
    email: str = Form(...),
    reference: str = Form("NM_007294.4"),
    alignment_type: str = Form("global"),
    variant_source: str = Form("clinvar"),
    population: str = Form(None),
    file: UploadFile = File(...)
):
    # Save uploaded file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    with open(temp_file.name, "wb") as f:
        shutil.copyfileobj(file.file, f)

    # Initialize and run analyzer
    analyzer = SequenceAnalyzer(reference=reference, email=email)
    if variant_source == "clinvar":
        analyzer.load_clinvar_data()

    result = analyzer.analyze_sample(temp_file.name)

    return JSONResponse(content=result)

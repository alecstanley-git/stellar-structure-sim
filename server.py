import subprocess
import json
from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

app = FastAPI()

app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/", response_class=HTMLResponse)
def read_root():
    with open("static/index.html") as f:
        return f.read()


@app.get("/simulate")
def simulate(mass: float = 1.0, x: float = 0.70, z: float = 0.02):
    y = max(0.0, 1.0 - x - z)
    # Run the rust cargo command
    cmd = [
        "cargo",
        "run",
        "--release",
        "--",
        "--mass",
        str(mass),
        "--x",
        str(x),
        "--y",
        str(y),
        "--z",
        str(z),
        "--json",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse output line by line, find the json part
        output = result.stdout.strip()
        lines = output.split("\n")
        # the json is the last line
        json_output = lines[-1]
        data = json.loads(json_output)
        return {"status": "success", "data": data}
    except subprocess.CalledProcessError as e:
        return {"status": "error", "message": e.stderr}


import uvicorn

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)

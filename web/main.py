from fastapi import FastAPI, File, UploadFile, Form, HTTPException, Depends, Cookie, Response, status
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse, FileResponse, RedirectResponse, JSONResponse
from fastapi.requests import Request
from fastapi.security import APIKeyCookie
import os
import sys
import tempfile
import logging
import shutil
from pathlib import Path
from typing import Optional

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add the parent directory to Python path so we can import mitocraft
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mitocraft
from mitocraft import MitocraftError, ReferenceBaseError, PipelineError, CommandError

# Create necessary directories
os.makedirs("web/static", exist_ok=True)
os.makedirs("final_output", exist_ok=True)

# Authentication settings
MITOEDIT_PASSWORD = os.getenv("MITOEDIT_PASSWORD")
if not MITOEDIT_PASSWORD:
    logger.error("MITOEDIT_PASSWORD environment variable is not set")
    raise ValueError("MITOEDIT_PASSWORD environment variable must be set to run the application")

CORRECT_PASSWORD = MITOEDIT_PASSWORD
COOKIE_NAME = "mitocraft_auth"
cookie_scheme = APIKeyCookie(name=COOKIE_NAME, auto_error=False)

logger.info("Using password from MITOEDIT_PASSWORD environment variable")

app = FastAPI(title="MitoEdit Web Interface")

# Mount static files and output directories
app.mount("/static", StaticFiles(directory="web/static"), name="static")
app.mount("/final_output", StaticFiles(directory="final_output"), name="final_output")

# Setup templates
templates = Jinja2Templates(directory="web/templates")

# Middleware to redirect HTTP to HTTPS when SSL is enabled
@app.middleware("http")
async def redirect_http_to_https(request: Request, call_next):
    # Check if SSL is enabled by looking for certificate files
    ssl_keyfile = os.getenv("SSL_KEYFILE", "certs/mitoedit.key")
    ssl_certfile = os.getenv("SSL_CERTFILE", "certs/mitoedit.pem")
    ssl_enabled = os.path.exists(ssl_keyfile) and os.path.exists(ssl_certfile)
    
    if ssl_enabled and request.url.scheme == "http":
        # Get the host and port from the request
        host = request.headers.get("host", "").split(":")[0]
        port = int(os.getenv("PORT", "443"))
        
        # Create the HTTPS URL
        url = request.url.replace(scheme="https", netloc=f"{host}:{port}")
        
        # Return a redirect response
        return RedirectResponse(url=str(url))
    
    # Continue with the request if not redirecting
    return await call_next(request)

# Authentication dependency
async def get_current_user(auth_token: Optional[str] = Cookie(None, alias=COOKIE_NAME)):
    if not auth_token or auth_token != CORRECT_PASSWORD:
        return None
    return {"authenticated": True}

# Login page
@app.get("/login", response_class=HTMLResponse)
async def login_page(request: Request):
    return templates.TemplateResponse(
        "login.html",
        {"request": request}
    )

# Login form submission
@app.post("/login")
async def login(response: Response, password: str = Form(...)):
    if password == CORRECT_PASSWORD:
        response = RedirectResponse(url="/", status_code=status.HTTP_302_FOUND)
        response.set_cookie(key=COOKIE_NAME, value=CORRECT_PASSWORD, httponly=True)
        return response
    else:
        return RedirectResponse(
            url="/login?error=Incorrect password",
            status_code=status.HTTP_302_FOUND
        )

# Main page - protected by authentication
@app.get("/", response_class=HTMLResponse)
async def home(request: Request, current_user = Depends(get_current_user)):
    if not current_user:
        return RedirectResponse(url="/login", status_code=status.HTTP_302_FOUND)
    
    return templates.TemplateResponse(
        "index.html",
        {"request": request}
    )

@app.post("/analyze")
async def analyze_sequence(
    request: Request,
    position: int = Form(...),
    reference_base: str = Form(...),
    mutant_base: str = Form(...),
    sequence_file: UploadFile = File(None),
    current_user = Depends(get_current_user)
):
    # Check authentication
    if not current_user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Authentication required"
        )
    
    # Initialize input_file
    input_file = None
    tmp_file = None

    try:
        # Get DNA sequence content
        if sequence_file and sequence_file.filename:
            logger.info(f"Processing uploaded sequence file: {sequence_file.filename}")
            suffix = '.fasta' if sequence_file.filename.endswith(('.fasta', '.fa')) else '.txt'
            content = await sequence_file.read()
            if isinstance(content, bytes):
                content = content.decode('utf-8')
            content = content.strip()
            logger.info(f"Uploaded file content length: {len(content)}")
            
            if not content:
                raise HTTPException(status_code=400, detail="Uploaded file is empty")
        else:
            logger.info("Using default mtDNA file")
            # Get content from default mtDNA file
            default_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'inputs', 'mito.txt')
            logger.info(f"Reading default mtDNA file: {default_file}")
            
            if not os.path.exists(default_file):
                raise HTTPException(status_code=500, detail="Default mtDNA file not found")
            
            with open(default_file, 'r') as src:
                content = src.read().strip()
                logger.info(f"Default mtDNA file content length: {len(content)}")

        if not content:
            raise HTTPException(status_code=500, detail="DNA sequence content is empty")

        # Create temporary file with the content
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.txt')
        tmp.write(content.encode('utf-8'))
        tmp.close()
        input_file = tmp.name
        tmp_file = input_file  # Store for cleanup
        logger.info(f"Created temporary file: {input_file}")

        # Verify the temporary file
        with open(input_file, 'r') as verify:
            verify_content = verify.read()
            logger.info(f"Verified temporary file content length: {len(verify_content)}")

        # Clean up and recreate output directories
        directories = ['running', 'final_output']
        for dir_name in directories:
            dir_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), dir_name)
            if os.path.exists(dir_path):
                logger.info(f"Cleaning up directory: {dir_path}")
                shutil.rmtree(dir_path)
            logger.info(f"Creating directory: {dir_path}")
            os.makedirs(dir_path, exist_ok=True)

        # Run mitocraft analysis
        try:
            logger.info(f"Starting analysis with position={position}, ref={reference_base}, mut={mutant_base}")
            # Override sys.argv with our parameters
            if sequence_file and sequence_file.filename:
                # If user uploaded a file, include the --input_file parameter
                args = [
                    'mitocraft.py',
                    '--input_file', input_file,
                    str(position),
                    reference_base.upper(),
                    mutant_base.upper()
                ]
            else:
                # If using default mtDNA file, don't include --input_file parameter
                args = [
                    'mitocraft.py',
                    str(position),
                    reference_base.upper(),
                    mutant_base.upper()
                ]
            logger.info(f"Running mitocraft with args: {args}")
            sys.argv = args
            mitocraft.main()
            logger.info("Analysis completed successfully")  
        except ReferenceBaseError as e:
            logger.error(f"Reference base error: {str(e)}")
            raise HTTPException(
                status_code=400, 
                detail=f"Reference base error: {str(e)}"
            )
        except PipelineError as e:
            logger.error(f"Pipeline error: {str(e)}")
            raise HTTPException(
                status_code=400, 
                detail=f"Pipeline error: {str(e)}"
            )
        except CommandError as e:
            logger.error(f"Command execution error: {str(e)}")
            raise HTTPException(
                status_code=500, 
                detail=f"Command execution error: {str(e)}"
            )
        except MitocraftError as e:
            logger.error(f"Mitocraft error: {str(e)}")
            raise HTTPException(
                status_code=500, 
                detail=f"Analysis error: {str(e)}"
            )
        except Exception as e:
            logger.error(f"Error during mitocraft analysis: {str(e)}")
            raise HTTPException(status_code=500, detail=f"Analysis error: {str(e)}")

        # Get the results file path
        results_file = f'final_output/final_{position}.xlsx'
        
        if os.path.exists(results_file):
            return {
                "status": "success",
                "results_file": results_file,
                "filename": f"final_{position}.xlsx"
            }
        else:
            raise HTTPException(status_code=500, detail="Analysis failed to generate results")

    except Exception as e:
        logger.error(f"Unexpected error: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Unexpected error: {str(e)}")
    finally:
        # Clean up temporary file if it was created
        if tmp_file and os.path.exists(tmp_file):
            try:
                os.unlink(tmp_file)
                logger.info(f"Cleaned up temporary file: {tmp_file}")
            except Exception as e:
                logger.error(f"Error cleaning up temporary file: {str(e)}")

if __name__ == "__main__":
    import uvicorn
    import os
    
    port = int(os.getenv("PORT", "8000"))
    
    # Check if SSL certificate and key files exist
    ssl_keyfile = os.getenv("SSL_KEYFILE", "certs/mitoedit.key")
    ssl_certfile = os.getenv("SSL_CERTFILE", "certs/mitoedit.pem")
    
    use_ssl = os.path.exists(ssl_keyfile) and os.path.exists(ssl_certfile)
    
    if use_ssl:
        logger.info(f"Starting HTTPS server on port {port} with SSL")
        uvicorn.run(
            app, 
            host="0.0.0.0", 
            port=port,
            ssl_keyfile=ssl_keyfile,
            ssl_certfile=ssl_certfile
        )
    else:
        logger.info(f"Starting HTTP server on port {port} without SSL")
        uvicorn.run(app, host="0.0.0.0", port=port)

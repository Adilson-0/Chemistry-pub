#importações das bibliotecas
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pesq import pesquisa
from pesq import Reacao

#criação do objeto FastAPI
app = FastAPI()

#configuração do CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

#rota para enviar a fórmula estrutural de uma molécula
@app.get('/molEstrutural/{molName}')
async def env_mol_e(molName):
    pesquisa(molName)
    return FileResponse('formula/Molestrutural.png')

#rota para enviar a fórmula bastão de uma molécula
@app.get('/molBastao')
async def env_mol_b():
    return FileResponse('formula/Molbastao.png')

#rota para enviar as informações e classificações de uma molécula
@app.get('/infos/mol/{molName}')
async def env_info_mol(molName):
    infos = pesquisa(molName)
    if(infos == "ErrorFind"):
        return {'Error':'Substância não encontrada!'}
    elif(infos == "ErrorFind_Draw"):
        return {'Error':'Não foi possível encontrar ou desenhar a substância!'}
    else:
        return infos
    
#rota para enviar a fórmula da reação química
@app.get('/reac/{str_reacao}')
async def env_reacao(str_reacao):
    Reacao(str_reacao)
    return FileResponse('formula/Reacao.png')

#rota para enviar as classificações da reação
@app.get('/infos/reac/{str_reacao}')
async def env_info_reac(str_reacao):
    infos = Reacao(str_reacao)
    if(len(infos) == 1):
        return {'Error':'Reação não encontrada!'}
    else:
        return {'class':infos[0], 'tipo':infos[1]}
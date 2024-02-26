from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions
from deep_translator import GoogleTranslator
from pymongo import MongoClient
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import pubchempy
import links

def Cadeia(m1): #Função de classificação de moléculas dado um objeto Mol
    Classificacao = []
    Insaturacao = len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]=,#[#6]"))) #Calcula Insaturações
    HeteroAtoms = len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c][!C!c][C,c]"))) #Calcula heteroátomos
    Terciary_Quartenary_Carbons = len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c][C,c]([C,c])[C,c]"))) + len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c]-[C,c]([C,c])([C,c])-[C,c]"))) #carbonos terciários e quaternários
    Aromatic_Rings = rdMolDescriptors.CalcNumAromaticRings(m1) #Número de núcleos aromáticos
    if (Aromatic_Rings > 0): # se tiver núcleo aromático
            Classificacao += ["Aromática"]
    elif (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c]@[C,c]"))) > 0): #se tiver carbonos dentro de uma cadeia fechada
        Classificacao += ["Alicíclica"]
    else:
        Classificacao += ["Acíclica"]
    if (Classificacao[0] == "Aromática"): #se for aromático
        Classificacao += ["Insaturada"] #já é insaturada
        if (Terciary_Quartenary_Carbons > 0): #se tiver carbono terciário ou quarternário
            Classificacao += ["Ramificada"]
        else:
            Classificacao += ["Normal"]
        if (HeteroAtoms > 0): #se tiver heteroátomo
            Classificacao += ["Heterogênea"]
        else:
            Classificacao += ["Homogênea"]
        if (Aromatic_Rings == 1): #se tiver um núcleo aromático
            Classificacao += ["Mononuclear"]
        elif (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c]:[C,c](:[C,c]):[C,c]"))) > 0): #se tiver mais de um núcleo e os carbonos dos núcleos estiverem ligados
            Classificacao += ["Condensada e Polinuclear"]
        else:
            Classificacao += ["Isolada e Polinuclear"]
    else: #caso da cadeia não aromática
        if (Insaturacao > 0): #se tiver insaturação
            Classificacao += ["Insaturada"]
        else:
            Classificacao += ["Saturada"]
        if (Classificacao[0] == "Alicíclica"): #se a cadeia for fechada
            if (Terciary_Quartenary_Carbons > 0): #se tiver carbonos terciários ou quarternários
                Classificacao += ["Ramificada"]
            else:
                Classificacao += ["Normal"]
            if (HeteroAtoms > 0): #se tiver heteroátomo
                Classificacao += ["Heterocíclica"]
            else:
                Classificacao += ["Homocíclica"]
        else: #caso cadeia não é fechada nem aromática
            if (Terciary_Quartenary_Carbons > 0): # se tiver carbonos terciários e quarternários sabendo que não é fechada nem aromática
                Classificacao += ["Ramificada"]
            else:
                Classificacao += ["Normal"]
            if (HeteroAtoms > 0): #se tiver heteroatomo
                Classificacao += ["Heterogênea"]
            else:
                Classificacao += ["Homogênea"]
    if (len(Classificacao) == 4): # caso não seja aromático faltará uma Classificacao(mononuclear,polinuclear e isolado,polinuclear e condensado), então:
        Classificacao += ["-"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c](=[O,o])[OX2H1,oX2H1]"))) > 0): #A partir daqui são as funções orgânicas
        Classificacao += ["Ácido Carboxílico"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[$([#6][CH1]=O),$([CH2]=O)]"))) > 0):
        Classificacao += ["Aldeído"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c][C,c](=[O,o])[C,c]"))) > 0):
        Classificacao += ["Cetona"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("C=[CX3]([OX2H,oX2H])"))) > 0):
        Classificacao += ["Enol"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("c:c([OX2H,oX2H]):c"))) > 0):
        Classificacao += ["Fenol"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[CX4][OX2H,oX2H]"))) > 0):
        Classificacao += ["Álcool"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[CX3H0,cX3H0](=[O,o])[OX2H0,oX2H0]")))-2*len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8]=[#6][#8][#6]=[#8]"))) > 0): #Este cálculo do éster tira todos os casos da estrutura de um anidrido (O-C=O=C-O) do cálculo do primeiro padrão (C(=O)O(R)R'), e como esse padrão do éster se repete duas vezes no anidrido, é tirado duas vezes.
        Classificacao += ["Éster"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c][O,o][C,c]")))-len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[CX3H0,cX3H0](=[O,o])[OX2H0,oX2H0]"))) > 0): #Do primeiro cálculo do padrão do éter (R-O-R') tira-se o padrão do éster (R-C(=O)O-R)
        Classificacao += ["Éter"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[N+,n+](=[O,o])[O-,o-]"))) > 0):
        Classificacao += ["Nitrocomposto"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[NX3,nX3]"))) - len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c](=[O,o])[NX3,nX3]"))) - len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[N+,n+](=[O,o])[O-,o-]"))) > 0): #Na conta do padrão da amina (R-NH2 ou R-N(R')H ou R-N(R')R'') tira-se o padrão da amida (R-C(=O)-NH2 ou C(=O)-N(R')H ou C(=O)-N(R')R'') e o do nitrocomposto (N+(=O)-O)
        Classificacao += ["Amina"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c](=[O,o])[NX3,nX3]")))  > 0):
        Classificacao += ["Amida"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c]~[F,Cl,Br,I]"))) > 0):
        Classificacao += ["Haleto"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c]#[N,n]"))) > 0):
        Classificacao += ["Cianeto"]
    if (len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8]=[#6][#8][#6]=[#8]"))) > 0):
        Classificacao += ["Anidrido"]
    if (len(m1.GetAtoms()) - len(m1.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]"))) == 0): #No hidrocarboneto são contados apenas os átomos que não são hidrogênios (pois o GetAtoms não captura hidrogênios a menos que a molécula tenha sido alterada pelo AddHs) e retira os carbonos.
        Classificacao += ["Hidrocarboneto"]
    if (len(Classificacao) > 6): #Caso tenha mais de uma função
        Classificacao = Classificacao[:5:] + ["Função Mista"]
    return Classificacao

def FormulaEstrutural(Mol): #Modifica os atomLabels do carbono para o grafo da molécula parecer uma fórmula estrutural condensada dado um objeto Mol
    for atom in Mol.GetAtoms(): #Cpatura cada átomo da molécula
        if (atom.GetTotalNumHs() > 1 and atom.GetSymbol() == "C"): #Caso número de H maior que 1 e o elemento tenha o símbolo do carbono.
            atom.SetProp("atomLabel","{}H{}".format(atom.GetSymbol(),atom.GetTotalNumHs())) #Então o atomLabel será CHn
        elif (atom.GetTotalNumHs() == 1 and atom.GetSymbol() == "C"): #Senão se o número de H igual a 1 e o símbolo do átomo é C
            atom.SetProp("atomLabel","{}H".format(atom.GetSymbol())) #Então o atomLabel é CH
        elif (atom.GetSymbol() == "C"): #Senão se o símbolo for C
            atom.SetProp("atomLabel",atom.GetSymbol()) #Então o atomLabel é C

def PesquisaPubchem(Nome): #Pesquisa no banco de dados do PubChem dado o nome da substância
    translator = GoogleTranslator(source="pt", target="en") #Objeto de tradução
    Traduzido = translator.translate(Nome + " (químico)") #O nome é traduzido, a substring " (químico)" está na string para evitar erros de tradução como propino -> bribe(propina)
    Traduzido = Traduzido[:len(Traduzido)-11:] #Particionamento para retirar a substring " (chemical)"
    try: #Tenta pela primeira vez pesquisar a substância
        comp = pubchempy.get_compounds(Traduzido,"name")[0] #Pesquisa composto pelo nome, retorna um vetor com os compostos(Compound object) pegos(apenas um)
        return comp.to_dict(properties=["isomeric_smiles"])["isomeric_smiles"] #Converte o objeto Compound em um dicionário com as propriedades definidas. Nessa caso, apenas o SMILES isomérico que é retornado.
    except: #Algumas substâncias no banco de dados do Pubchem não têm o último hífen. Exemplo: 5-ethyl-3-methyl-4-propyl-octane está no PubChem como 5-ethyl-3-methyl-4-propyloctane
        Traduzido = Traduzido[::-1].replace("-","",1)[::-1] #Inverte a string e retira o primeiro(que era o último) e desinverte-a 
        comp = pubchempy.get_compounds(Traduzido,"name")[0] #Faz a busca como na linha anterior
        return comp.to_dict(properties=["isomeric_smiles"])["isomeric_smiles"] #Retorna o Smiles Isomérico

def pesquisa(mol): #Pesquisa a substância pelo nome no banco de dados primário
    Smiles = ""
    Molecule = ""
    m = ""
    Hidrocarbonet = mol.lower().strip() #Baixar a caixa da letra e tirar os espaços no final da string
    try:
        Cluster = MongoClient(links.LINK_DB) #Cria objeto MongoClient conectando-se ao Cluster
        Database = Cluster.get_database("chemistry") #Obtém o banco de dados como objeto
        Collection = Database.get_collection("substancias") #Obtém a collections substancias como objeto
        Smiles = Collection.find({"nome" : Hidrocarbonet},{"smiles" : 1})[0]["smiles"] #Faz uma operação de consulta no banco de dados por uma substância com certo nome e retorna o SMILES dela
        m = rdmolfiles.MolFromSmiles(Smiles) #Cria um objeto Mol a partir do SMILES
        Molecule = Chem.rdMolDescriptors.CalcMolFormula(m) #Calcula a fórmula molecular
        Info = {1 : m, "Nome do composto" : Hidrocarbonet, "Fórmula Molecular" : Molecule, "Classificacoes" : Cadeia(m)} 
        DesenharMol(m)
        return Info
    except: #Caso dê errado
        try: #Processo similar, porém com o PubChem
            Smiles = PesquisaPubchem(Hidrocarbonet)
            m = rdmolfiles.MolFromSmiles(Smiles)
            Molecule = Chem.rdMolDescriptors.CalcMolFormula(m)
            Info = {1 : m, "Nome do composto" : Hidrocarbonet, "Fórmula Molecular" : Molecule, "Classificacoes" : Cadeia(m)}
            DesenharMol(m)
            return Info
        except:
            Error = 'ErrorFind'
            return Error

def DesenharMol(m1): #Cria a imagem da molécula
    if (rdmolfiles.MolToSmiles(m1) != "C"): #Identifica se o SMILES da substância é o mesmo do metano
        Draw.MolToFile(m1,links.DIRETORIO_FORMULA_BASTAO,size=(500,400)) #Cria a fórmula bastão(naturalmente feita pelo comando) da substância com dimensões 500x400
        FormulaEstrutural(m1) #Modifica os atomLabels dos carbonos para eles aparecerem no grafo da molécula
        drawer = rdMolDraw2D.MolDraw2DCairo(500,400) #Cria um objeto MolDraw2DCairo
        drawer.drawOptions().fixedFontSize = 20 #Chama as opções do objeto e configura o tamanho fixo da fonte para 20 pixels
        drawer.DrawMolecule(m1) #Renderiza a imagem da molécula
        drawer.FinishDrawing() #Adiciona os últimos bits
        drawer.WriteDrawingText(links.DIRETORIO_FORMULA_ESTRUTURAL) #Guarda a imagem no local especificado
    else: #Se for
        m1 = rdmolops.AddHs(m1) #Adiciona os hidrogênios ao grafo da molécula
        I = Image.new("RGB",(500,500),color=(255,255,255)) #Cria um objeto image no modo RGB(pixels de 3xi bit) com a cor branca.
        Im = ImageDraw.Draw(I) #Cria o objeto para desenhar na imagem criada na linha anterior
        Font = ImageFont.truetype("fonte\ARIAL.TTF",size=25) #Define a fonte como Arial e o tamanho das letras como 25 pixels
        Im.text((25,250),"Como o metano tem apenas um carbono \nnão é possível fazer sua fórmula bastão.", font=Font,fill = (0,0,0)) #Escreve a string nas coordenadas (25, 250) e com cor preta
        I.save(links.DIRETORIO_FORMULA_BASTAO) #Guarda a imagem no local especificado
        for atom in m1.GetAtoms(): #Captura os átomos da molécula
            atom.SetProp("atomLabel",atom.GetSymbol()) #Configura o atomLabel para ficar igual ao símbolo do átomo
        Draw.MolToFile(m1,links.DIRETORIO_FORMULA_ESTRUTURAL,size=(500,400)) #Guarda a imagem no local especificado
        

def ClassificaReacao(Reacao): #Classifica uma reação dado um objeto ChemicalReaction
    Reagentes = []
    Produtos = []
    for Reactants in Reacao.GetReactants(): #Obtém os reagentes da reação
        K = rdmolops.SanitizeMol(Reactants) #Verifica as valências, calcula outras representações para a molécula(kekulization), configura aromaticidade e conjugação.
        T = True
        Result = Cadeia(Reactants) #Classifica o reagente
        if ("Saturada" in Result and "Hidrocarboneto" in Result and "Acíclica" in Result): # A partir daqui usa as classificações do reagente e outros padrões adicionais para definí-lo
            Reagentes += ["Alcano"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[$([#17][#17]),$([#35][#35]),$([#9][#9]),$([#53][#53])]"))) > 0):
            Reagentes += ["Halogênio"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#7+](=[#8])([#8H1])[#8-]"))) > 0):
            Reagentes += ["Ácido nítrico"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8H1][#16](=[#8])(=[#8])[#8H1]"))) > 0):
            Reagentes += ["Ácido sulfúrico"]
        if ("Aromática" in Result and "Hidrocarboneto" in Result):
            Reagentes += ["Hidrocarboneto aromático"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[$([#17H1]),$([#35H1]),$([#9H1]),$([#53H1])]"))) > 0):
            Reagentes += ["Haleto de hidrogênio"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8]=[#6][!#6;!#1]"))) > 0):
            Reagentes += ["Grupo Acila"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8H2]"))) > 0):
            Reagentes += ["Água"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]=[#6]"))) == 1 and len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]#[#6]"))) == 0 and "Hidrocarboneto" in Result):
            Reagentes += ["Alceno"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#1][#1]"))) > 0):
            Reagentes += ["Hidrogênio"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]#[#6]"))) + len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]=[#6]"))) >= 1 and "Hidrocarboneto" in Result):
            Reagentes += ["Hidrocarboneto Insaturado"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]@[C]"))) > 0 and "Hidrocarboneto" in Result and "Saturada" in Result):
            Reagentes += ["Ciclano"]
        if ("Álcool" in Result):
            Reagentes += ["Álcool"]
        if ("Ácido Carboxílico" in Result):
            Reagentes += ["Ácido Carboxílico"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8-][#8+]=[#8]"))) > 0):
            Reagentes += ["Ozônio"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8X0]"))) > 0):
            Reagentes += ["Oxigênio atômico"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#19+].[#8H1-]"))) > 0):
            Reagentes += ["Hidróxido de potássio"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6H3][#8H1]"))) > 0):
            Reagentes += ["Metanol"]
        if (len(Reactants.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#1X0]"))) > 0):
            Reagentes += ["Hidrogênio atômico"]
        if ("Aldeído" in Result):
            Reagentes += ["Aldeído"]
        if ("Cetona" in Result):
            Reagentes += ["Cetona"]
        if ("Amida" in Result):
            Reagentes += ["Amida"]
        if ("Nitrocomposto" in Result):
            Reagentes += ["Nitrocomposto"]
        if ("Haleto" in Result):
            Reagentes += ["Haleto"]

    for products in Reacao.GetProducts(): #Obtém os produtos de uma reação
        K = rdmolops.SanitizeMol(products) #O mesmo que ocorre com cada reagente ocorre com cada produto
        Result = Cadeia(products) #Classificação do produto
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#1][#1]"))) > 0): #Mesma lógica dos reagentes
            Produtos += ["Hidrogênio"]
        if ("Insaturada" in Result and "Hidrocarboneto" in Result):
            Produtos += ["Composto insaturado"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[$([#17H1]),$([#35H1]),$([#9H1]),$([#53H1])]"))) > 0):
            Produtos += ["Haleto de hidrogênio"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8H2]"))) > 0):
            Produtos += ["Água"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6][#16](=[#8])(=[#8])[#8H1]"))) > 0):
            Produtos += ["Ácido sulfônico"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("c[$([#17]),$([#9]),$([#35]),$([#53])]"))) > 0):
            Produtos += ["Haleto Aromático"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("c[#7+](=[#8])[#8-]"))) > 0):
            Produtos += ["Nitrocomposto Aromático"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("c[#16](=[#8])(=[#8])[#8]"))) > 0):
            Produtos += ["Ácido sulfônico Aromático"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("c[#6]"))) > 0 and "Hidrocarboneto" in Result):
            Produtos += ["Hidrocarboneto Aromático Ramificado"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("c[#6](=[#8])[#6]"))) > 0):
            Produtos += ["Cetona Aromática"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]=[C]"))) == 1 and "Hidrocarboneto" in Result):
            Produtos += ["Alceno"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[$([#17][#17]),$([#35][#35]),$([#9][#9]),$([#53][#53])]"))) > 0):
            Produtos += ["Halogênio"]
        if ("Hidrocarboneto" in Result and "Saturada" in Result):
            Produtos += ["Alcano"]
        if ("Haleto" in Result and "Saturada" in Result):
            Produtos += ["Haleto Saturado"]
        if ("Haleto" in Result and "Insaturada" in Result and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]=[C]"))) == 1):
            Produtos += ["Haleto Insaturado"]
        if ("Nitrocomposto" in Result and not "Nitrocomposto Aromático" in Produtos):
            Produtos += ["Nitrocomposto"]
        if ("Álcool" in Result and "Saturada" in Result):
            Produtos += ["Álcool Saturado"]
        if ("Haleto" in Result and "Insaturada" in Result and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]=[C]"))) > 1):
            Produtos += ["Haleto Insaturado com mais de uma insaturação"]
        if ("Álcool" in Result and "Insaturada" in Result and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]=[C]"))) > 1):
            Produtos += ["Álcool Insaturado com mais de uma insaturação"]
        if ("Álcool" in Result and "Insaturada" in Result and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]=[C]"))) == 1):
            Produtos += ["Álcool Insaturado"]
        if ("Haleto" in Result and "Saturada" in Result and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]@[C][$([#17]),$([#9]),$([#35]),$([#53])]"))) > 0):
            Produtos += ["Haleto cíclico e saturado"]
        if ("Haleto" in Result and "Acíclica" in Result):
            Produtos += ["Haleto aberto"]
        if ("Haleto" in Result and "Alicíclica" in Result and "Insaturada" in Result):
            Produtos += ["Haleto Insaturado com mais de uma insaturação"]
        if ("Anidrido" in Result):
            Produtos += ["Anidrido"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#19].[$([#17]),$([#9]),$([#53]),$([#35])]"))) > 0):
            Produtos += ["Hidróxido com haleto"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#19][$([#17]),$([#35]),$([#9]),$([#53])]"))) > 0):
            Produtos += ["Haleto de potássio"]
        if ("Aldeído" in Result):
            Produtos += ["Aldeído"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8H1][#8H1]"))) > 0):
            Produtos += ["Peróxido de hidrogênio"]
        if ("Cetona" in Result and not ("Cetona Aromática" in Result)):
            Produtos += ["Cetona"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8H1][#6][#6][#8H1]"))) > 0):
            Produtos += ["Diol vicinal"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]#[#6]"))) + len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#6]=[#6]"))) >= 1 and "Hidrocarboneto" in Result):
            Produtos += ["Hidrocarboneto Insaturado"]
        if ("Saturada" in Result and "Hidrocarboneto" in Result):
            Produtos += ["Hidrocarboneto Saturado"]
        if ("Ácido Carboxílico" in Result):
            Produtos += ["Ácido Carboxílico"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[#8]=[#6]=[#8]"))) > 0):
            Produtos += ["Dióxido de Carbono"]
        if ("Éter" in Result):
            Produtos += ["Éter"]
        if ("Amina" in Result):
            Produtos += ["Amina"]
        if ("Hidrocarboneto" in Result):
            Produtos += ["Hidrocarboneto"]
        if ("Álcool" in Result):
            Produtos += ["Álcool"]
        if ("Cianeto" in Result):
            Produtos += ["Cianeto"]
        if ("Função Mista" in Result and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c](=[O,o])[OX2H1,oX2H1]"))) > 0 and len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C,c][C,c](=[O,o])[C,c]"))) > 0):
            Produtos += ["Ácido Carboxílico","Cetona"]
        if (len(products.GetSubstructMatches(rdmolfiles.MolFromSmarts("[C]@[C]"))) > 0 and "Hidrocarboneto" in Result and "Saturada" in Result):
            Produtos += ["Ciclano"]
    Classificacao = []
#A partir daqui as classificações e definições dos produtos e reagentes são usados para definir a classificação da reação
    if (("Alcano" in Reagentes or "Hidrocarboneto aromático" in Reagentes) and "Halogênio" in Reagentes and ("Haleto Saturado" in Produtos or "Haleto Aromático" in Produtos) and "Haleto de hidrogênio" in Produtos):
        Classificacao += ["Substituição","Halogenação"]
    if (("Alcano" in Reagentes or "Hidrocarboneto aromático" in Reagentes) and "Ácido nítrico" in Reagentes and ("Nitrocomposto" in Produtos or "Nitrocomposto Aromático" in Produtos) and "Água" in Produtos):
        Classificacao += ["Substituição","Nitração"]
    if (("Alcano" in Reagentes or "Hidrocarboneto aromático" in Reagentes) and "Ácido sulfúrico" in Reagentes and ("Ácido sulfônico" in Produtos or "Ácido sulfônico Aromático" in Produtos) and "Água" in Produtos):
        Classificacao += ["Substituição","Sulfonação"]
    if ("Hidrocarboneto aromático" in Reagentes and "Haleto" in Reagentes and "Haleto de hidrogênio" in Produtos and "Hidrocarboneto Aromático Ramificado" in Produtos):
        Classificacao += ["Substituição","Alquilação"]
    if ("Hidrocarboneto aromático" in Reagentes and "Grupo Acila" in Reagentes and "Haleto de hidrogênio" in Produtos and "Cetona Aromática" in Produtos):
        Classificacao += ["Substituição","Alcilação"]
    if ("Haleto" in Reagentes and "Água" in Reagentes and ("Álcool Saturado" in Produtos or "Álcool Insaturado com mais de uma insaturação" in Produtos or "Álcool insaturado" in Produtos) and "Haleto de hidrogênio" in Produtos):
        Classificacao += ["Substituição","Hidrólise Alcalina"]
    if ("Alceno" in Reagentes and "Hidrogênio" in Reagentes and "Alcano" in Produtos):
        Classificacao += ["Adição","Hidrogenação"]
    if ("Alceno" in Reagentes and "Halogênio" in Reagentes and "Haleto Saturado" in Produtos):
        Classificacao += ["Adição","Halogenação"]
    if ("Alceno" in Reagentes and "Haleto de hidrogênio" in Reagentes and "Haleto Saturado" in Produtos):
        Classificacao += ["Adição","Hidro-halogenação"]
    if ("Alceno" in Reagentes and "Água" in Reagentes and "Álcool Saturado" in Produtos):
        Classificacao += ["Adição","Hidratação"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Hidrogênio" in Reagentes and "Hidrocarboneto Insaturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Hidrogenação parcial"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Hidrogênio" in Reagentes and "Hidrocarboneto Saturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Hidrogenação total"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Haleto de hidrogênio" in Reagentes and "Haleto Insaturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Hidro-Halogenação parcial"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Haleto de hidrogênio" in Reagentes and "Haleto Saturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Hidro-Halogenação total"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Halogênio" in Reagentes and "Haleto Insaturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Halogenação parcial"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Halogênio" in Reagentes and "Haleto Saturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Halogenação total"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Água" in Reagentes and "Álcool Insaturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Hidratação parcial"]
    if ("Hidrocarboneto Insaturado" in Reagentes and "Água" in Reagentes and "Álcool Saturado" in Produtos and len(Classificacao) == 0):
        Classificacao += ["Adição","Hidratação total"]
    if ("Hidrocarboneto aromático" in Reagentes and "Halogênio" in Reagentes and "Haleto cíclico e saturado" in Produtos):
        Classificacao += ["Adição","Halogenação"]
    if ("Hidrocarboneto aromático" in Reagentes and "Hidrogênio" in Reagentes and "Ciclano" in Produtos):
        Classificacao += ["Adição","Hidrogenação"]
    if ("Ciclano" in Reagentes and "Hidrogênio" in Reagentes and "Alcano" in Produtos):
        Classificacao += ["Adição","Hidrogenação"]
    if ("Ciclano" in Reagentes and "Halogênio" in Reagentes and ("Haleto aberto" in Produtos or ("Haleto cíclico e saturado" in Produtos and "Haleto de hidrogênio" in Produtos))):
        Classificacao += ["Adição","Halogenação"]
    if (("Alcano" in Reagentes or "Álcool" in Reagentes or "Amina" in Reagentes) and "Hidrogênio" in Produtos and ("Composto insaturado" in Produtos or "Cianeto" in Produtos or "Cetona" in Produtos or "Aldeído" in Produtos)):
        Classificacao += ["Eliminação","Desidrogenação"]
    if ("Haleto" in Reagentes and "Alceno" in Produtos and not "Haleto de hidrogênio" in Produtos):
        Classificacao += ["Eliminação","De-halogenação"]
    if ("Haleto" in Reagentes and "Alceno" in Produtos and "Haleto de hidrogênio" in Produtos):
        Classificacao += ["Eliminação","Desidrohalogenação"]
    if ("Álcool" in Reagentes and "Alceno" in Produtos and "Água" in Produtos):
        Classificacao += ["Eliminação","Desidratação intramolecular"]
    if ("Álcool" in Reagentes and "Éter" in Produtos and "Água" in Produtos):
        Classificacao += ["Eliminação","Desidratação intermolecular"]
    if (Reagentes.count("Ácido Carboxílico") > 1 and "Anidrido" in Produtos and "Água" in Produtos):
        Classificacao += ["Eliminação","Desidratação intermolecular"]
    if (Reagentes.count("Ácido Carboxílico") == 1 and "Anidrido" in Produtos and "Água" in Produtos):
        Classificacao += ["Eliminação","Desidratação intramolecular"]
    if ("Haleto" in Reagentes and "Hidróxido de potássio" in Reagentes and "Hidróxido com haleto" in Produtos and "Alceno" in Produtos and "Água" in Produtos):
        Classificacao += ["Eliminação","-"]
    if ("Alceno" in Reagentes and "Ozônio" in Reagentes and "Água" in Reagentes and ("Aldeído" in Produtos or "Cetona" in Produtos) and "Peróxido de hidrogênio" in Produtos):
        Classificacao += ["Oxidação","Ozonólise"]
    if ("Alceno" in Reagentes and "Oxigênio atômico" in Reagentes and ("Diol vicinal" in Produtos or "Aldeído" in Produtos or "Cetona" in Produtos)):
        Classificacao += ["Oxidação","Oxidação Branda"]
    if (("Alceno" in Reagentes or "Ciclano" in Reagentes) and "Oxigênio atômico" in Reagentes and ("Ácido Carboxílico" in Produtos or "Cetona" in Produtos)):
        Classificacao += ["Oxidação","Oxidação Energética"]
    if ("Álcool" in Reagentes and "Oxigênio atômico" in Reagentes and "Cetona" in Produtos and "Água" in Produtos):
        Classificacao += ["Oxidação","Oxidação Energética"]
    if ("Álcool" in Reagentes and "Oxigênio atômico" in Reagentes and "Aldeído" in Produtos and "Água" in Produtos):
        Classificacao += ["Oxidação","Oxidação Parcial"]
    if ("Álcool" in Reagentes and "Oxigênio atômico" in Reagentes and "Ácido Carboxílico" in Produtos and "Água" in Produtos):
        Classificacao += ["Oxidação","Oxidação Total"]
    if ("Álcool" in Reagentes and "Oxigênio atômico" in Reagentes and "Dióxido de Carbono" in Produtos and "Água" in Produtos):
        Classificacao += ["Oxidação","-"]
    if ("Álcool" in Reagentes and "Haleto de hidrogênio" in Reagentes and "Hidrocarboneto" in Produtos and "Halogênio" in Produtos and "Água" in Produtos):
        Classificacao += ["Redução","-"]
    if ("Ácido Carboxílico" in Reagentes and "Hidrogênio atômico" in Reagentes and "Aldeído" in Produtos and "Água" in Produtos):
        Classificacao += ["Redução","Redução Parcial"]
    if ("Ácido Carboxílico" in Reagentes and "Hidrogênio atômico" in Reagentes and "Álcool" in Produtos):
        Classificacao += ["Redução","Redução Total"]
    if ("Aldeído" in Reagentes and "Hidrogênio atômico" in Reagentes and "Álcool" in Produtos):
        Classificacao += ["Redução","-"]
    if ("Cetona" in Reagentes and "Hidrogênio atômico" in Reagentes and "Álcool" in Produtos):
        Classificacao += ["Redução","-"]
    if ("Amida" in Reagentes and "Hidrogênio atômico" in Reagentes and "Amina" in Produtos and "Água" in Produtos):
        Classificacao += ["Redução","-"]
    if ("Nitrocomposto" in Reagentes and "Hidrogênio atômico" in Reagentes and "Amina" in Produtos and "Água" in Produtos):
        Classificacao += ["Redução","-"]

    return Classificacao

def Reacao(R): #Pesquisa de reação química
    R = E = R.lower().strip() #Baixa a caixa da letra da reação química
    Cluster = MongoClient(links.LINK_DB)
    Database = Cluster.get_database("chemistry")
    reacoes = Database.get_collection("reacoes") #Collection reacoes como objeto collection
    substancias = Database.get_collection("substancias")
    Outras = Database.get_collection("outras") #Collection outras como objeto collection
    Rm = R.split() #Separação da string pelos espaços (caso o primeiro bloco dê erro a solução seria montar a reação a partir dessa string com os nomes dos reagentes,catalisadores ou solventes e produtos separados)
    try: 
        R = reacoes.find({"reação" : R.replace("_"," ")},{"smiles" : 1})[0]["smiles"] #Procura a reação no base de dados do MongoDB e retorna o SMILES dela
        Reacao = rdChemReactions.ReactionFromSmarts(R, useSmiles=True) #Cria um objeto ChemicalReaction(a reação) a partir do SMILES
        Reactants = Reacao.GetReactants()
        Agents = Reacao.GetAgents() 
        Products = Reacao.GetProducts()
        Reaction = rdChemReactions.ChemicalReaction() #Criar uma ChemicalReaction vazia
        for Reactant in Reactants:
            K = rdmolops.SanitizeMol(Reactant) 
            FormulaEstrutural(Reactant) #Para deixar os atomLabels dos carbonos visíveis na reação
            Reaction.AddReactantTemplate(Reactant) #Adiciona o reagente na ChemicalReaction vazia
        for Agent in Agents:
            K = rdmolops.SanitizeMol(Agent) #Mesma função sendo aplicada aos catalisadores ou produtos
            FormulaEstrutural(Agent) #Para deixar os atomLabels dos carbonos visíveis na reação
            Reaction.AddAgentTemplate(Agent) #Adiciona o catalisador ou solvente na ChemicalReaction vazia
        for Product in Products:
            K = rdmolops.SanitizeMol(Product)
            FormulaEstrutural(Product) #Para deixar os atomLabels dos carbonos visíveis na reação
            Reaction.AddProductTemplate(Product) #Adiciona o produto na ChemicalReaction vazia
        drawer = rdMolDraw2D.MolDraw2DCairo(2000,800)
        drawer.DrawReaction(Reaction) #Renderiza a imagem da reação
        drawer.FinishDrawing()
        drawer.WriteDrawingText(links.DIRETORIO_REACAO)
        return ClassificaReacao(Reaction) #Classifica a reação
    except: #Montagem da reação por partes
        ReacaoR = rdChemReactions.ChemicalReaction()
        i = 0
        T = True
        P = False
        while i < len(Rm) and T:
            if (Rm[i]._contains_(">") and Rm[i] != "->" and Rm[i] != "+"): #Caso a substância esteja no espaço dos catalisadores ou solventes e não seja um + nem ->
                Catal = Rm[i] #Catalisadores
                Catal = Catal[1:len(Catal)-1:] #tira o '-' e o ">"
                CatalList = Catal.split("+") #Separa as substâncias pelo '+'
                for var2 in range(len(CatalList)): 
                    try: #Caso N não exista só será adicionada 1 vez (mesmo processo)
                        CatalList[var2] = Outras.find({"nome" : CatalList[var2].replace("_"," ")},{"smiles" : 1})[0]["smiles"]
                        Agent = rdmolfiles.MolFromSmiles(CatalList[var2])
                        FormulaEstrutural(Agent)
                        ReacaoR.AddAgentTemplate(Agent)
                    except: #Caso o catalisador ou solvente não exista no banco de dados do MongoDB
                        try: #Faz o mesmo processo mas com o PubChem
                            CatalList[var2] = PesquisaPubchem(CatalList[var2])
                            Agent = rdmolfiles.MolFromSmiles(CatalList[var2])
                            FormulaEstrutural(Agent)
                            ReacaoR.AddAgentTemplate(Agent)
                        except: #Caso dê errado a variável T vai ser Falsa
                            T = False
                P = True #No fim do for o código chegará aos produtos
            elif P == False and Rm[i] != "+" and Rm[i] != "->": #Caso a substância esteja no espaço de reagentes (pois P é falso) e não seja um + nem ->
                try: #Adiciona apenas 1 vez o produto (mesmo processo do catalisador ou solvente, mas dessa vez a pesquisa ocorre em substancias)
                    v = substancias.find({"nome" : Rm[i].replace("_"," ")},{"smiles" : 1})[0]["smiles"]
                    Reactant = rdmolfiles.MolFromSmiles(v)
                    FormulaEstrutural(Reactant)
                    ReacaoR.AddReactantTemplate(Reactant)
                except: #Caso dê errado
                    try: #Mesmo processo mas na collection Outras
                        v = Outras.find({"nome" : Rm[i].replace("_"," ")},{"smiles" : 1})[0]["smiles"]
                        Reactant = rdmolfiles.MolFromSmiles(v)
                        FormulaEstrutural(Reactant)
                        ReacaoR.AddReactantTemplate(Reactant)
                    except: #Caso dê errado
                        try: #PubChem (mesmo processo do catalisador ou solkvente)
                            v = PesquisaPubchem(Rm[i].replace("_"," "))
                            Reactant = rdmolfiles.MolFromSmiles(v)
                            FormulaEstrutural(Reactant)
                            ReacaoR.AddReactantTemplate(Reactant)
                        except: #Caso dê errado T é Falso
                            T = False
            elif Rm[i] != "+" and Rm[i] != "->": #Caso a substância esteja no espaço dos produtos (pois P é True), não seja + nem ->
                try: #Mesmo processo dos produtos
                    v = substancias.find({"nome" : Rm[i].replace("_"," ")},{"smiles" : 1})[0]["smiles"]
                    Product = rdmolfiles.MolFromSmiles(v)
                    FormulaEstrutural(Product)
                    ReacaoR.AddProductTemplate(Product)
                except: #Caso dê errado
                    try: #Mesmo caso dos produtos
                        v = Outras.find({"nome" : Rm[i].replace("_"," ")},{"smiles" : 1})[0]["smiles"]
                        Product = rdmolfiles.MolFromSmiles(v)
                        FormulaEstrutural(Product)
                        ReacaoR.AddProductTemplate(Product)
                    except: #Caso dê errado
                        try: #Pesquisa no PubChem
                            v = PesquisaPubchem(Rm[i].replace("_",""))
                            Reactant = rdmolfiles.MolFromSmiles(v)
                            FormulaEstrutural(Product)
                            ReacaoR.AddProductTemplate(Product)
                        except: #Caso dê errado
                            T = False #T é falso
            elif (Rm[i] == "->"): #Caso não haja catalisador ou solvente apenas a seta existirá e quando ela passar...
                P = True #Chega-se aos produtos
            i += 1
        if T == True: #Se T for True é porque todos os reagentes, catalisadores ou solventes, produtos foram encontrados e adicionados
            drawer = rdMolDraw2D.MolDraw2DCairo(2000,800) #Mesmo processo que começa na linha 432
            drawer.DrawReaction(ReacaoR)
            drawer.FinishDrawing()
            drawer.WriteDrawingText(links.DIRETORIO_REACAO)
            return ClassificaReacao(ReacaoR)
        else:
            return ["Reação não encontrada."]
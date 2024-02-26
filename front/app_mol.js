//div com a imagem da fórmula estrutural
const div_molEst = document.getElementById('formulaE');
//div com a imagem da fórmula bastão
const div_molBas = document.getElementById('formulaB');
//div com as informações da molécula
const div_infos = document.getElementById('inf_resp');
//botão de pesquisa
const button_pesq_mol = document.getElementById('button_mol');
//caixa de texto que recebe o nome da molécula
const input_pesq = document.getElementById('inp_comp');
//fórmulário de busca por molécula
const form = document.getElementById('pesq_form');

//evento que será ativado caso clique no botão
button_pesq_mol.addEventListener('click', () => {
    //rota que envia a fórmula estrutural
    var url_molEst = 'http://127.0.0.1:8000/molEstrutural/'+input_pesq.value;
    //rota que envia as informações
    var url_infos = 'http://127.0.0.1:8000/infos/mol/'+input_pesq.value;
    //rota que envia a fórmula bastão
    var url_molBas = 'http://127.0.0.1:8000/molBastao/';

    //requisição para a rota da fórmula estrutural
    fetch(url_molEst, {method:'GET'})
    .then(resp => {
        //transforma o arquivo recebido em um Blob
        return resp.blob();
    })
    .then(img => {
        //cria uma URL para o arquivo recebido
        const img_url = URL.createObjectURL(img);   
        div_molEst.innerHTML = `<img id="img_molE" src="${img_url}" />`;
    })
    .catch(error => console.log(error)); //caso ocorra um erro...

    //requisição para a rota da fórmula bastão
    fetch(url_molBas, {method:'GET'})
    .then(resp => {
        //transforma o arquivo recebido em um Blob
        return resp.blob();
    })
    .then(img => {
        //cria uma URL para o arquivo recebido
        const img_url = URL.createObjectURL(img);
        div_molBas.innerHTML = `<img id="img_molB" src="${img_url}" />`;
    })
    .catch(error => console.log(error)); //caso ocorra um erro...

    //requisição para a rota das informações
    fetch(url_infos, {method:'GET'})
    .then(resp => {
        //transforma os dados recebidos em um JSON
        return resp.json()
    })
    .then(infos => {
        if(Object.keys(infos).length == 4){
            //coleta dos dados do JSON
            var molName = infos['Nome do composto'];
            var formulaMol = infos['Fórmula Molecular'];
            var classi = infos['Classificacoes'];

            //elemento HTML que será inserido na div
            var informs = `<ul>
                <li><strong>nome do composto:</strong> ${molName}</li>
                <li><strong>fórmula molecular:</strong> ${formulaMol}</li>
                <li><strong>classificações:</strong>
                    <ul style="margin-left: 5%;">
                        <div style="display: flex;">
                            <div>
                                <li>${classi[0]}</li>
                                <li>${classi[1]}</li>
                                <li>${classi[2]}</li>
                            </div>
                            <div style="margin-left: 20%;">
                                <li>${classi[3]}</li>
                                <li>${classi[4]}</li>
                                <li>${classi[5]}</li>
                            </div>
                        </div>
                    </ul>
                </li>
            </ul>`;

            div_infos.innerHTML = informs;
        }
        
        //caso a substância não tenha sido encontrada
        else{
            var erro = infos['Error'];
            var erro_msg = `<p style="color: rgb(255, 74, 74); font-size: 120%;">Erro: ${erro}</p>`;
            div_infos.innerHTML = erro_msg;
        }
    })
    .catch(error => console.log(error)); //caso ocorra algum erro...});
});


//evento que será ativado caso pressione enter no input
form.addEventListener('submit', (evento) => {
    evento.preventDefault();
});


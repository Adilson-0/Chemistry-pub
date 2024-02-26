//botão de adição de reagente
const button_add_reag = document.getElementById("button_reac_reag");
//botão de adição de catalisador/solvente
const button_add_catal = document.getElementById("button_reac_catal");
//botão de adição de produto
const button_add_prod = document.getElementById("button_reac_prod");
//botão que inicia o processo de pesquisa
const button_pesq = document.getElementById("button_reac");
//botão para cancelamento de pesquisa
const button_cancel_pesq = document.getElementById("button_cancel_reac");
//caixa de texto para os reagentes
const input_reag = document.getElementById("inp_reagente");
//caixa de texto para os catalisadores/solventes 
const input_catal = document.getElementById("inp_catalizador");
//caixa de texto para os produtos
const input_prod = document.getElementById("inp_produto");
//div que contém a imagem da reação
const div_reac = document.getElementById("reac");
//div que apresenta as classificações da reação
const div_infos_reac = document.getElementById("inf_resp_reac");
//div que exibe os reagentes adicionados
const div_reags = document.getElementById("reags");
//div que exibe os reagentes adicionados
const div_catals = document.getElementById("catals");
//div que exibe os reagentes adicionados
const div_prods = document.getElementById("prods");

//variável da parte dos reagentes
var string_reag = "";
//variável da parte dos catalisadores/solventes
var string_catal = "";
//variável da parte dos produtos
var string_prod = "";
//variável da string de busca (completa)
var string_pesq = "";

//número de adições de reagentes
var num_reag = 0;
//número de adições de catalisadores/solventes
var num_catal = 0;
//número de adições de produtos
var num_prod = 0;

button_add_reag.addEventListener('click', ()=>{
    if(num_reag >= 1){ 
        string_reag += " + " + input_reag.value.trim().replace(/ /g, "_");       
        num_reag++;
        div_reags.innerHTML += `<p>${input_reag.value}</p>`;
        input_reag.value = "";
    }else{
        string_reag += input_reag.value.trim().replace(/ /g, "_");
        num_reag++;
        div_reags.innerHTML += `<p>${input_reag.value}</p>`;
        input_reag.value = "";
    }
});

button_add_catal.addEventListener('click', ()=>{
    if(num_catal >= 1){
        string_catal += "+" + input_catal.value.trim().replace(/ /g, "_");
        num_catal++;
        div_catals.innerHTML += `<p>${input_catal.value}</p>`;
        input_catal.value = "";
    }else{
        string_catal += input_catal.value.trim().replace(/ /g, "_");
        num_catal++;
        div_catals.innerHTML += `<p>${input_catal.value}</p>`;
        input_catal.value = "";
    }
});

button_add_prod.addEventListener('click', ()=>{
    if(num_prod >= 1){
        string_prod += " + " + input_prod.value.trim().replace(/ /g, "_");
        num_prod++;
        div_prods.innerHTML += `<p>${input_prod.value}</p>`;
        input_prod.value = "";
    }else{
        string_prod += input_prod.value.trim().replace(/ /g, "_");
        num_prod++;
        div_prods.innerHTML += `<p>${input_prod.value}</p>`;
        input_prod.value = "";
    }
});

button_cancel_pesq.addEventListener('click', ()=>{
    string_catal = "";
    string_pesq  = "";
    string_prod  = "";
    string_reag  = "";
    num_catal = 0;
    num_prod  = 0;
    num_reag  = 0;

    div_catals.innerHTML = `<p></p>`
    div_reags.innerHTML = `<p></p>`;
    div_prods.innerHTML = `<p></p>`;
});

button_pesq.addEventListener('click', ()=>{
    string_pesq = "";
    //string de busca completa
    string_pesq = string_reag + " -" + string_catal + "> " + string_prod;

    //rotas para as requisições
    var url_reac_img = 'http://127.0.0.1:8000/reac/'+string_pesq;
    var url_infos_reac = 'http://127.0.0.1:8000/infos/reac/'+string_pesq;

    //requisição de imagem da reação
    fetch(url_reac_img, {method:'GET'})
    .then(resp => {
        //transforma o arquivo recebido em um Blob
        return resp.blob();
    })
    .then(img => {
        //cria uma URL para o arquivo recebido
        const img_url = URL.createObjectURL(img);
        div_reac.innerHTML = `<img id="img_reac" src="${img_url}" />`;
    })
    .catch(error => console.log(error));//caso ocorra um erro...

    //requisição das classificações da reação
    fetch(url_infos_reac, {method:'GET'})
    .then(resp => {
        //transforma os dados recebidos em um JSON
        return resp.json()
    })
    .then(infos => {
        if(Object.keys(infos).length == 2){
            var classi = infos['class'];
            var tipo = infos['tipo'];

            //elemento HTML que será inserido na div
            var informs = `<ul>
            <li><strong>Classificação da reação:</strong> ${classi}</li>
            <li><strong>Tipo de reação:</strong> ${tipo}</li>
        </ul>`;

            div_infos_reac.innerHTML = informs;
        }
        //caso a reação não tenha sido encontrada
        else{
            var erro = infos['Error'];
            var erro_msg = `<p style="color: rgb(255, 74, 74); font-size: 120%;">Erro: ${erro}</p>`;
            div_infos_reac.innerHTML = erro_msg;
            div_reac.innerHTML = `<img src="" alt="erro"/>`;
        }
    })
    .catch(error => console.log(error));//caso ocorra um erro...
    
    //reset das variáveis de contagem e das partes da sring de busca
    string_catal = "";
    string_pesq = "";
    string_prod = "";
    string_reag = "";
    url_infos_reac = "";
    url_reac_img = "";
    num_catal = 0;
    num_prod = 0;
    num_reag = 0;
    div_catals.innerHTML = `<p></p>`;
    div_reags.innerHTML = `<p></p>`;
    div_prods.innerHTML = `<p></p>`;
})


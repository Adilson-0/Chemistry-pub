<div align="center">
  <h1>Chemistry</h1>
  <p>Aplicação web para realizar pesquisas por fórmulas de moléculas e reações orgânicas. Projeto desenvolvido para TCC do curso técnico de informática do IFPB - Campina Grande</p>
</div>

**Alunos:**
<ul>
  <li>Adilson fernandes da Silva Filho</li>
  <li>Gleryston Lucas da Silva Marques</li>
</ul>

<hr>

<h2>Como funciona o site</h2>
<p>Para pesquisar uma molécula orgânica basta ir até a parte "pesquisa de molécula", digitar o nome do composto na caixa de texto, e clicar no botão "iniciar pesquisa". Caso a molécula for encontrada, será exibida a sua fórmula estrutural, fórmula bastão, e algumas informações adicionais. Caso ocorra algum erro, ou a molécula não for encontrada, será exibido uma mensagem de erro na área "informações sobre a substância".</p> 
<div align="center">
  <img src="https://docs.google.com/uc?id=1YayOJ6OIjuhM28EOgyDf2BYlnrEz8YUq" height="400" width="600"/>
</div>
<br>
<p>Já para pesquisar por reações químicas, basta ir até a parte "pesquisa de reação" e escrever os reagentes, solventes/catalisadores e produtos em seus respectivos campos. Cada vez que um composto for adicionado, deve clicar no botão "adicionar" para integrá-lo à reação que será pesquisada, após isso o campo será limpo, permitindo a inserção de outros compostos.<br>Quando todas substâncias forem adicionadas, clique em "iniciar pesquisa". Se a reação não for encontrada, será exibido uma mensagem de erro, semelhante a da seção anterior.</p>
<div align="center">
  <img src="https://docs.google.com/uc?id=1tYGiO2WwJ_6AifVi1ewIiBxjW1MZnVyv" height="550" width="600"/>
</div>

<h2>Como executar a API</h2>
<p>Para executar a API é necessário instalar as bibliotecas utilizadas nos códigos python e em seguida executar os seguintes comandos no terminal:</p>
<br>
<ul>
  <li>cd api</li>
  <li>uvicorn main:app --reload</li>
</ul>

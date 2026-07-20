# viralfetch — Especificação de Projeto

> Documento de contexto para o Claude Code. Leia por inteiro antes de escrever código.
> Implemente **por fases**, na ordem indicada. Não pule para fases posteriores.
> Implemente tudo em inglês.

---

## 1. O que é

`viralfetch` é uma ferramenta de linha de comando em Python para consultar e baixar
taxonomia viral, metadados e sequências, combinando três fontes:

| Fonte | Conteúdo | Acesso |
|---|---|---|
| **VMR** (ICTV Virus Metadata Resource) | hierarquia taxonômica + isolados exemplares + acessos GenBank/RefSeq | **local**, TSV embarcado |
| **NCBI E-utilities** | metadados de sequências, sequências (nt/aa), taxonomia NCBI | **remoto**, sob demanda |
| **ICTV Report** (`ictv.global/report/chapter/...`) | texto descritivo dos capítulos por família/ordem | **remoto**, sob demanda (scraping) |

Princípio central: **o VMR é o índice local; tudo o mais é buscado sob demanda e cacheado.**

---

## 2. Estado atual do repositório

```
├── data
│   └── VMR_MSL41.v1.20260320.tsv
├── LICENSE
└── README.md
```

O arquivo TSV é a release MSL41 do VMR. Ele é a fonte de verdade para taxonomia.
**Nunca modifique esse arquivo.** Ele é read-only; qualquer normalização gera artefatos derivados.

---

## 3. Restrições técnicas inegociáveis

Estas regras existem por motivos de política de uso e de robustez. Não as contorne.

### NCBI

- Use **apenas E-utilities** (`eutils.ncbi.nlm.nih.gov/entrez/eutils/*.fcgi`).
  **Nunca** faça scraping de páginas HTML do NCBI (`ncbi.nlm.nih.gov/nuccore/...`).
- Todo request deve incluir os parâmetros `tool=viralfetch` e `email=<do usuário>`.
  O e-mail vem de `$NCBI_EMAIL` ou de arquivo de config. **Se não houver e-mail
  configurado, o comando falha com mensagem explicativa** — não invente um valor padrão.
- Chave de API opcional, lida de `$NCBI_API_KEY`.
- Rate limit: **3 req/s sem chave, 10 req/s com chave.** Implemente um limitador central
  no cliente, não `sleep()` espalhado pelo código.
- Use **POST**, não GET, em todos os endpoints (listas de acessos estouram o limite de URL).
- Lotes de no máximo **200 acessos** por request.
- Falha parcial é normal: se você pede 200 acessos e voltam 197, **reporte os 3 ausentes**
  ao usuário. Não falhe silenciosamente nem estoure exceção.
- Faça argumentos --store-ncbi-email e --store-ncbi-apikey para salvar os dados do usuário, caso ele queira.

### ICTV

- `User-Agent`: `viralfetch/<versão> (+<url do repo>; <email>)`
- Delay mínimo de **1 segundo** entre requisições.
- Respeite `robots.txt` (verifique em `https://ictv.global/robots.txt`).
- Conteúdo é **CC BY 4.0** — toda saída de texto deve incluir a atribuição
  (autores do capítulo + citação/DOI quando disponível).

### Cache

- Diretório: `platformdirs.user_cache_dir("viralfetch")`
- Sequências e metadados de acessos: **cache permanente** (acessos são imutáveis).
- Textos de capítulos: **TTL de 30 dias**.
- Comandos `viralfetch cache info` e `viralfetch cache clear [--texts|--seqs]`.
- Flag global `--no-cache` para forçar refetch.

---

## 4. Arquitetura

```
viralfetch/
├── data/
│   └── VMR_MSL41.v1.20260320.tsv        # read-only
├── src/viralfetch/
│   ├── __init__.py
│   ├── cli.py               # Typer: só define comandos e chama serviços
│   ├── config.py            # e-mail, api key, paths, resolução de config
│   ├── models.py            # dataclasses: Taxon, Isolate, Sequence, Chapter, ...
│   ├── vmr.py               # carga e indexação do TSV
│   ├── accessions.py        # parser do campo de acesso (texto livre → lista)
│   ├── ncbi.py              # cliente E-utilities (esummary/efetch/elink/esearch)
│   ├── ictv.py              # fetch + parse dos capítulos → markdown
│   ├── cache.py             # cache em disco com TTL
│   └── render/
│       ├── __init__.py      # despacha entre rich e json
│       ├── rich_.py         # todas as tabelas/painéis Rich
│       └── json_.py         # serialização JSON
├── tests/
│   ├── fixtures/            # HTML e respostas do NCBI congelados
│   └── ...
├── pyproject.toml
├── LICENSE
└── README.md
```

### Regra arquitetural mais importante

**Nenhuma função de lógica de negócio imprime nada.** Todas retornam
dataclasses ou dicts. A camada `render/` é a única que escreve no stdout.

Isso é o que faz a flag `--json` funcionar sem duplicar código:

```python
# cli.py
resultado = servico.consultar_taxon(nome)     # retorna dados puros
render.emitir(resultado, formato=ctx.formato) # rich ou json
```

Se você se pegar chamando `console.print()` dentro de `vmr.py` ou `ncbi.py`, está errado.

### Dependências

```
typer          # CLI
rich           # renderização
requests       # HTTP
platformdirs   # paths de cache/config
beautifulsoup4 # parsing do ICTV
lxml           # parser rápido para o bs4
```

Nada além disso sem justificativa. **Não use pandas** — o VMR cabe em memória
como lista de dataclasses e pandas é peso morto num CLI (tempo de import).

---

## 5. Superfície de comandos

Flags globais: `--json`, `--no-cache`, `--verbose`, `--email`, `--api-key`.

### 5.1 Taxonomia pura (VMR local)

```bash
viralfetch tax <nome>
viralfetch tax Coronaviridae
viralfetch tax "Betacoronavirus pandemicum"
```

Mostra a linhagem completa do táxon (realm → espécie), o rank detectado,
e um resumo dos isolados associados se for espécie.

Busca deve ser **case-insensitive** e aceitar match parcial com sugestão
(`Você quis dizer: ...`) quando não houver match exato.

### 5.2 Taxonomia comparada (ICTV vs NCBI)

```bash
viralfetch tax <nome> --compare-ncbi
```

Fluxo:
1. Obtém a linhagem ICTV do VMR local.
2. Pega um acesso representativo da espécie.
3. `elink` de `nuccore` → `taxonomy` para obter o **taxid do NCBI**.
4. `efetch` em `db=taxonomy` para obter a linhagem NCBI.
5. Renderiza as duas linhagens **lado a lado**, destacando divergências.

Divergências são esperadas e comuns (o NCBI atrasa em relação ao ICTV).
Trate isso como o produto do comando, não como erro.

### 5.3 Membros de um táxon

```bash
viralfetch members <taxon> [--rank <rank>] [--count]
viralfetch members Coronaviridae --rank genus
viralfetch members Riboviria --rank family --count
```

Lista os táxons filhos em qualquer nível abaixo do informado.
`--count` mostra só as quantidades agregadas por rank.
Puramente local, sem rede.

### 5.4 Dados NCBI de uma espécie

```bash
viralfetch seq <espécie> --meta
viralfetch seq <espécie> --fasta [-o arquivo.fa]
viralfetch seq <espécie> --gb
```

- `--meta` → `esummary` em `db=nuccore`. Rápido, ~1 KB por acesso.
  Campos a exibir: `accessionversion`, `organism`, `slen`, `moltype`, `biomol`,
  `topology`, `completeness`, `sourcedb`, `updatedate`.
- `--fasta` → `efetch` com `rettype=fasta`.
- `--gb` → `efetch` com `rettype=gb` (registro completo com features).

Sem flag de formato, o padrão é `--meta`.

### 5.5 Dados NCBI de um táxon inteiro

```bash
viralfetch seq --taxon <taxon> --meta          # resumo agregado
viralfetch seq --taxon <taxon> --fasta -o out.fa
viralfetch seq --taxon Filoviridae --meta
```

Com `--meta` num táxon, mostre um **agregado**: quantas espécies, quantos
isolados, quantos acessos, quebra por `moltype`, e quantos são RefSeq.
Isso é o que permite ao usuário decidir se vale baixar.

**Sempre confirme antes de downloads grandes.** Se o táxon tem mais de 500
acessos, mostre o total e peça confirmação (a menos que `--yes` seja passado).

### 5.6 Filtro por tipo de molécula

```bash
viralfetch seq --taxon <taxon> --moltype ssRNA --fasta
viralfetch seq <espécie> --biomol mRNA
viralfetch seq <espécie> --protein --fasta
```

**Atenção — este é o ponto de maior confusão do projeto:**

- `moltype` e `biomol` (genomic, mRNA, cRNA, ...) são campos de `db=nuccore`.
  Filtre-os **sobre o resultado do `esummary`**, localmente.
- **Aminoácidos NÃO são um filtro de `nuccore`.** Proteínas vivem em `db=protein`.
  Para obtê-las: `elink` de `nuccore` → `protein`, depois `efetch` em `db=protein`.
  Implemente isso como caminho separado, acionado por `--protein`.

Documente essa distinção no `--help` do comando, porque não é óbvia para o usuário.

### 5.7 Texto do ICTV Report

```bash
viralfetch text <família>
viralfetch text Coronaviridae [--section summary] [--raw]
```

- Busca `https://ictv.global/report/chapter/<slug>/<slug>`.
- Converte o HTML do conteúdo principal em Markdown.
- Renderiza com `rich.markdown.Markdown`.
- `--raw` emite o Markdown puro (para redirecionar a arquivo).
- Sempre inclua no rodapé: autores, citação/DOI, e a nota de licença CC BY 4.0.

**Cuidados no parser:**
- Preserve o **itálico dos nomes científicos** — é obrigatório em nomenclatura viral.
- O menu de navegação global se repete em toda página. Restrinja os seletores ao
  container de conteúdo principal; nunca use `soup.find_all("a")` sem escopo.
- Preserve tabelas de características (converta para tabela Markdown).
- Falhe **alto e explicitamente** se um campo esperado sumir. Não grave `None`
  em silêncio — isso mascara quebras causadas por mudança de tema no site.

### 5.8 Utilitários

```bash
viralfetch update            # verifica se há VMR mais novo em ictv.global/vmr
viralfetch cache info
viralfetch cache clear [--texts|--seqs]
viralfetch config            # mostra e-mail, api key (mascarada), paths
```

---

## 6. Parsing de acessos (crítico)

A coluna de acesso do VMR é **texto livre**, não identificador limpo:

```
NC_045512                              → simples
RNA1: NC_003615; RNA2: NC_003616       → segmentado
DNA-A: X15656; DNA-B: X15657           → geminivírus
AB012345 (partial)                     → com anotação
```

Vírus segmentados (Reoviridae, Orthomyxoviridae, Nanoviridae) tornam o
tratamento de segmentos obrigatório. Ignorá-lo perde genomas inteiros.

O parser deve retornar `list[tuple[segmento | None, acesso]]` e normalizar
o modelo de dados para **uma linha por acesso**, não por espécie.

Ao carregar o VMR, **conte e exponha quantas linhas produziram zero acessos.**
Esse número é o indicador de qualidade do parser. Inclua um comando ou flag
de diagnóstico que liste esses casos.

---

## 7. Convenções de saída

- **Rich por padrão** em tudo: tabelas para listas, painéis para detalhes,
  árvores (`rich.tree`) para linhagens, barras de progresso para downloads.
- `--json` produz JSON puro no stdout, **sem nenhuma decoração**, pronto para `jq`.
  Nesse modo, mensagens de progresso e avisos vão para **stderr**.
- Erros sempre em stderr, exit code diferente de zero.
- Se stdout não for um TTY, desabilite cores automaticamente (o Rich já faz isso,
  mas confirme que o comportamento está correto em pipes).

---

## 8. Fases de implementação

Implemente e valide cada fase antes de passar à seguinte.

### Fase 1 — Fundação
- `pyproject.toml`, estrutura de diretórios, entry point `viralfetch`.
- `models.py` com as dataclasses.
- `vmr.py`: carrega o TSV, constrói índices por nome e por rank.
- `accessions.py` com o parser e testes.
- **Aceite:** `python -c "from viralfetch.vmr import load; print(len(load().species))"` funciona.

### Fase 2 — Taxonomia local
- Comandos `tax` (sem `--compare-ncbi`) e `members`.
- Camada `render/` com Rich e JSON.
- **Aceite:** `viralfetch tax Coronaviridae` e `viralfetch members Coronaviridae --rank genus`
  funcionam nos dois formatos de saída.

### Fase 3 — Cliente NCBI
- `ncbi.py` com rate limiting, retry com backoff exponencial, POST, lotes.
- `cache.py`.
- Comandos `seq` com `--meta`, `--fasta`, `--gb`.
- **Aceite:** `viralfetch seq "Betacoronavirus pandemicum" --meta` retorna dados reais.

### Fase 4 — NCBI avançado
- `seq --taxon` com agregação e confirmação de downloads grandes.
- Filtros `--moltype` / `--biomol`.
- Caminho de proteína via `elink` (`--protein`).
- `tax --compare-ncbi`.

### Fase 5 — Textos do ICTV
- `ictv.py`: fetch, parse, conversão para Markdown.
- Comando `text`.
- Fixtures de teste com 5 capítulos congelados (incluindo `Coronaviridae`
  e um geminivírus, que têm estruturas diferentes).

### Fase 6 — Acabamento
- `update`, `cache`, `config`.
- Autocompletar de shell para nomes de táxon (Typer suporta nativamente).
- README com exemplos reais de uso.

---

## 9. Testes

- **Nenhum teste deve fazer requisição de rede.** Use fixtures congeladas em
  `tests/fixtures/` (respostas JSON do NCBI, HTML do ICTV).
- Teste o parser de acessos contra os casos difíceis listados na seção 6.
- Teste de regressão do parser de capítulos: se o HTML do ICTV mudar,
  o teste deve quebrar de forma óbvia.
- Marque testes de integração com `@pytest.mark.network` e exclua do padrão.

---

## 10. Não faça

- Não use pandas.
- Não faça scraping de HTML do NCBI.
- Não faça nada em português (tudo em inglês);
- Não imprima de dentro da lógica de negócio.
- Não invente um e-mail padrão para os parâmetros do NCBI.
- Não modifique o TSV do VMR.
- Não engula exceções de rede — reporte falhas parciais ao usuário.
- Não implemente cache "inteligente" com invalidação complexa. TTL simples basta.
- Não adicione TUI, servidor web, ou qualquer coisa fora do escopo acima.

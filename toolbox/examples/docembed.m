function [E,info]=docembed(docs,varargin)
%docembed computes neural embeddings of text documents, as unit vectors, calling Python
%
%<a href="matlab: docsearchFS('docembed')">Link to the help function</a>
%
%   An embedding is a vector of p numbers (p from a few hundred to a few
%   thousand) which a neural network, trained on a very large collection of
%   texts, associates to a document: documents which talk about the same
%   subject have vectors which point in similar directions. The similarity
%   of two documents is the cosine of the angle between their vectors, and
%   the length of the vectors carries no information. For this reason the
%   rows of the output are normalised to length 1 and are points of the
%   unit hypersphere, the space of tclustvmf. MATLAB has no function to
%   compute this type of embedding for generic text: docembed writes the
%   documents in a file, runs a Python script (which is included at the end
%   of this file, so no other file is needed) and reads the embeddings back.
%
%  Required input arguments:
%
%        docs: documents. Cell array of character vectors (or string array)
%              of length n. Each element is the text of one document.
%
%  Optional input arguments:
%
%     backend: type of embedding. Character.
%              'sbert' sentence transformer (default). A transformer
%                      network which reads the text (about the first 256 to
%                      512 word pieces) and returns one vector for the whole
%                      text. Needs: pip install sentence-transformers
%              'glove' average of GloVe word vectors, weighted with idf
%                      (static word embeddings, 300 dimensions). Needs:
%                      pip install gensim
%              'spacy' average of the word vectors of a spaCy model. Needs:
%                      pip install spacy and a model (en_core_web_md).
%              'auto'  'sbert' if it is available, otherwise 'glove'.
%                 Example - 'backend','glove'
%                 Data Types - char
%       model: model of the backend. Character. For 'sbert' the default is
%              'all-MiniLM-L6-v2' (384 dimensions, 90 MB, fast). Other
%              choices are 'all-MiniLM-L12-v2' (384), 'all-mpnet-base-v2'
%              (768, the best of this family), 'paraphrase-MiniLM-L6-v2'
%              (384), 'BAAI/bge-small-en-v1.5' (384), 'BAAI/bge-base-en-v1.5'
%              (768), 'intfloat/e5-small-v2' (384),
%              'paraphrase-multilingual-MiniLM-L12-v2' (384, 50 languages),
%              or the path of a folder with a model saved on disk (offline
%              machines). The first time a model is downloaded from
%              huggingface.co and stored in the cache of the user
%              (~/.cache/huggingface, or the folder in the environment
%              variable HF_HOME). For 'glove' the default is
%              'glove-wiki-gigaword-300' (others for gensim: 'glove-wiki-
%              gigaword-100', 'fasttext-wiki-news-subwords-300'); for
%              'spacy' it is 'en_core_web_md'.
%                 Example - 'model','all-mpnet-base-v2'
%                 Data Types - char
%   pythonexe: Python to use. Character. Path of the Python executable. The
%              default is the one selected by pyenv.
%                 Example - 'pythonexe','/usr/bin/python3'
%                 Data Types - char
%   cachefile: file in which the embeddings are saved. Character. If the
%              file exists and was produced for the same documents, backend
%              and model, the embeddings are read from it and Python is not
%              run. The default is '' (no cache). With many documents the
%              computation takes minutes or hours: use it.
%                 Example - 'cachefile','docembeddings.mat'
%                 Data Types - char
%         msg: display information on the screen. Scalar. If msg=1
%              (default) the Python command and its messages are displayed.
%                 Example - 'msg',0
%                 Data Types - double
%
%  Output:
%
%           E: embeddings. Matrix of size n-by-p, of class double. Each row
%              has length 1.
%        info: information on the computation. Structure with fields
%              backend (the backend actually used), model, pythonexe,
%              seconds (time of the Python run) and fromcache (true if the
%              embeddings were read from cachefile).
%
% More About:
%
% HOW MATLAB CALLS PYTHON. The "MATLAB Engine API for Python" is the
% interface in the opposite direction (Python calls MATLAB). To call Python
% from MATLAB the tools are pyenv (selects the installation of Python),
% py.module.function (calls a function) and pyrunfile (runs a script).
% docembed runs the script in a separate process with system, using the
% executable of pyenv: large libraries such as torch do not have to be loaded
% in the MATLAB process, where they sometimes clash with the libraries of
% MATLAB. The script, which can also be run from a shell, is the text
% between the lines %%%PYTHON-BEGIN and %%%PYTHON-END at the end of this file.
%
% PRECISION. The embeddings are returned in double precision. Neural
% networks work in single precision and the .mat file written by Python
% would be of class single; the functions of the suite (for example
% tclustvmf) require double.
%
% LONG TEXTS. Transformer models read only the first part of the text
% (256 to 512 word pieces, according to the model). For help pages this is
% the title and the description of the function. For long documents split
% the text in chunks, embed the chunks and average the normalised vectors.
%
% See also: doccorpus, tclustvmf, sphereprojplot
%
% References:
%
%   Reimers N., Gurevych I. (2019), Sentence-BERT: Sentence embeddings
%   using Siamese BERT-networks, Proceedings of EMNLP-IJCNLP 2019,
%   pp. 3982-3992.
%
%   Mu J., Viswanath P. (2018), All-but-the-top: simple and effective
%   postprocessing for word representations, ICLR 2018.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('docembed')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Embeddings of the help pages of FSDA and MATLAB, with the default
    % sentence transformer. Python with sentence-transformers is needed.
    [docs,docname,theme,themename]=doccorpus;
    E=docembed(docs,'cachefile','docembeddings.mat');
    disp(size(E))
%}
%
%{
    % The same with the average of GloVe word vectors (needs gensim).
    [docs,docname,theme,themename]=doccorpus;
    E=docembed(docs,'backend','glove');
%}

%% Beginning of code

options=struct('backend','sbert','model','','pythonexe','','cachefile','','msg',1);
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    if length(varargin)~=2*length(UserOptions)
        error('FSDA:docembed:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        error('FSDA:docembed:NonExistInputOpt','In total %d non-existent user options found.',length(WrongOptions));
    end
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
backend=options.backend;
model=options.model;
if isempty(model) && any(strcmp(backend,{'sbert','auto'}))
    model='all-MiniLM-L6-v2';
end
msg=options.msg;
cachefile=options.cachefile;

docs=cellstr(docs);
docs=docs(:);
n=numel(docs);
% a signature of the documents, to know if a cache file is still valid
sig=[n sum(cellfun(@(x) sum(double(x)),docs))];
requested=sprintf('%s|%s',backend,model);

info=struct('backend',backend,'model',model,'pythonexe','','seconds',0,'fromcache',false);

%% Cache
if ~isempty(cachefile) && isfile(cachefile)
    S=load(cachefile);
    if isfield(S,'sig') && isequal(S.sig,sig) && strcmp(S.requested,requested)
        E=double(S.E);
        info=S.info;
        info.fromcache=true;
        if msg
            fprintf('Embeddings read from %s (%s, %s)\n',cachefile,info.backend,info.model)
        end
        return
    end
end

%% Python executable
pythonexe=options.pythonexe;
if isempty(pythonexe)
    pe=pyenv;
    pythonexe=char(pe.Executable);
    if isempty(pythonexe)
        pythonexe='python3';
        warning('FSDA:docembed:NoPyenv', ...
            'pyenv does not list a Python: python3 is searched in the path of the system');
    end
end
info.pythonexe=pythonexe;

%% Write the documents and the Python script
tmpdir=tempname;
mkdir(tmpdir);
cleanup=onCleanup(@() rmdir(tmpdir,'s'));
jsonfile=fullfile(tmpdir,'docs.json');
matfile=fullfile(tmpdir,'emb.mat');
pyfile=fullfile(tmpdir,'docembed.py');
fid=fopen(jsonfile,'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(docs));
fclose(fid);
src=fileread([mfilename('fullpath') '.m']);
ib=strfind(src,'%%%PYTHON-BEGIN');
ie=strfind(src,'%%%PYTHON-END');
lines=regexp(src(ib(end):ie(end)-1),'\r?\n','split');
lines=lines(2:end);
fid=fopen(pyfile,'w','n','UTF-8');
for i=1:numel(lines)
    % each line of Python starts with the characters %: in this file
    if numel(lines{i})>=3
        fprintf(fid,'%s\n',lines{i}(4:end));
    else
        fprintf(fid,'\n');
    end
end
fclose(fid);

%% Run Python
cmd=sprintf('"%s" "%s" "%s" "%s" "%s" "%s" 2>&1',pythonexe,pyfile,jsonfile,matfile,backend,model);
if msg
    fprintf('Running: %s\n',cmd)
end
tic
[status,out]=system(cmd);
info.seconds=toc;
if msg
    fprintf('%s\n',out)
end
if status~=0
    error('FSDA:docembed:PythonFailed', ...
        ['Python failed. Message:\n%s\nThe packages required by the backend ' ...
        'must be installed in the Python which is used (see the help).'],out)
end

%% Read the embeddings
R=load(matfile);
E=double(R.E);
info.backend=char(R.backend);
info.model=char(R.model);
if msg
    fprintf('%d-by-%d embeddings (%s, %s) in %.0f seconds\n',size(E,1),size(E,2), ...
        info.backend,info.model,info.seconds)
end
if ~isempty(cachefile)
    save(cachefile,'E','sig','requested','info')
end

end

%%%PYTHON-BEGIN
%: #!/usr/bin/env python3
%: """Neural embeddings of a list of documents (script written and run by docembed.m).
%:
%: Usage: python docembed.py docs.json embeddings.mat backend model
%:   docs.json: JSON array of strings; backend: auto, sbert, glove or spacy.
%: Output .mat: E (n-by-p, unit-norm rows, float64), backend, model.
%: """
%: import sys, json, re
%: import numpy as np
%: from scipy.io import savemat
%:
%: STOP = set("""a about above after again all also am an and any are as at be because been before being below
%: between both but by can could did do does doing down during each few for from further had has have having he her
%: here hers him his how i if in into is it its itself just me more most my no nor not of off on once only or other
%: our out over own same she should so some such than that the their them then there these they this those through
%: to too under until up use used using very was we were what when where which while who whom why will with would
%: you your""".split())
%:
%: def tokens(t):
%:     return [w for w in re.findall(r"[A-Za-z]{3,}", t.lower()) if w not in STOP]
%:
%: def idf_weights(toks):
%:     n = len(toks); df = {}
%:     for t in toks:
%:         for w in set(t): df[w] = df.get(w, 0) + 1
%:     return {w: np.log(n / c) + 1 for w, c in df.items()}
%:
%: def emb_sbert(docs, model):
%:     from sentence_transformers import SentenceTransformer
%:     name = model or "all-MiniLM-L6-v2"
%:     m = SentenceTransformer(name)
%:     return np.asarray(m.encode(docs, batch_size=64, show_progress_bar=False)), name
%:
%: def emb_glove(docs, model):
%:     import gensim.downloader as api
%:     name = model or "glove-wiki-gigaword-300"
%:     kv = api.load(name)
%:     toks = [tokens(d) for d in docs]
%:     idf = idf_weights(toks)
%:     E = np.zeros((len(docs), kv.vector_size))
%:     for i, t in enumerate(toks):
%:         ws = [(w, idf[w]) for w in t if w in kv.key_to_index]
%:         if ws:
%:             E[i] = sum(c * kv[w] for w, c in ws) / sum(c for _, c in ws)
%:     return E, name
%:
%: def emb_spacy(docs, model):
%:     import spacy
%:     name = model or "en_core_web_md"
%:     nlp = spacy.load(name, disable=["parser", "ner"])
%:     toks = [tokens(d) for d in docs]
%:     idf = idf_weights(toks)
%:     E = np.zeros((len(docs), nlp.vocab.vectors_length))
%:     for i, t in enumerate(toks):
%:         ws = [(w, idf[w]) for w in t if nlp.vocab.has_vector(w)]
%:         if ws:
%:             E[i] = sum(c * nlp.vocab[w].vector for w, c in ws) / sum(c for _, c in ws)
%:     return E, name
%:
%: def embed(docs, backend="auto", model=""):
%:     if backend in ("auto", "sbert"):
%:         try:
%:             E, name = emb_sbert(docs, model); backend = "sbert"
%:         except Exception as e:
%:             if backend == "sbert": raise
%:             print("sbert not available (%s); falling back to glove" % type(e).__name__, file=sys.stderr)
%:             backend = "glove"; model = ""
%:     if backend == "glove": E, name = emb_glove(docs, model)
%:     elif backend == "spacy": E, name = emb_spacy(docs, model)
%:     elif backend != "sbert": raise ValueError("unknown backend " + backend)
%:     E = np.asarray(E, dtype=np.float64)       # double: required by the MATLAB functions
%:     E = E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-12)   # unit norm
%:     return E, backend, name
%:
%: if __name__ == "__main__":
%:     docs = json.load(open(sys.argv[1], encoding="utf-8"))
%:     if isinstance(docs, str): docs = [docs]      # MATLAB jsonencode of 1 document
%:     backend = sys.argv[3] if len(sys.argv) > 3 else "auto"
%:     model = sys.argv[4] if len(sys.argv) > 4 else ""
%:     E, backend, name = embed(docs, backend, model)
%:     savemat(sys.argv[2], {"E": E, "backend": backend, "model": name})
%:     print("saved %d x %d embeddings (%s, %s)" % (E.shape[0], E.shape[1], backend, name))
%%%PYTHON-END

%FScategory:CLUS-Sphere

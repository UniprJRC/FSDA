// Configurazione
const EXTERNAL_URL = 'https://rosa.unipr.it/FSDA/index.html'; // Sostituisci con l'URL reale
const TARGET_DIV_ID = 'content_container'; // ID del div da cui vuoi estrarre il contenuto

// Funzione per caricare contenuto esterno
async function loadExternalContent() {
    const contentContainer = document.getElementById('content-container');
    
    try {
        // Nota: per problemi CORS, potresti dover usare un proxy o configurare il server esterno
        // Opzione 1: Usa un proxy CORS (per testing)
        const proxyUrl = `https://api.allorigins.win/get?url=${encodeURIComponent(EXTERNAL_URL)}`;
        
        const response = await fetch(proxyUrl);
        console.log(response);
        const data = await response.json();
        
        // Crea un elemento temporaneo per parsare l'HTML
        const tempDiv = document.createElement('div');
        tempDiv.innerHTML = data.contents;
        
        // Trova il div specifico
        const targetContent = tempDiv.querySelector(`#${TARGET_DIV_ID}`);
        
        if (targetContent) {
            contentContainer.innerHTML = targetContent.innerHTML;
        } else {
            contentContainer.innerHTML = '<div class="error">Contenuto non trovato nella pagina esterna.</div>';
        }
        
    } catch (error) {
        console.error('Errore durante il caricamento:', error);
        contentContainer.innerHTML = `
            <div class="error">
                <h3>Errore durante il caricamento del contenuto</h3>
                <p>Non è stato possibile caricare il contenuto dalla pagina esterna.</p>
                <p>Dettagli: ${error.message}</p>
            </div>
        `;
    }
}

// Gestione click menu laterale
document.addEventListener('DOMContentLoaded', function() {
    // Carica il contenuto iniziale
    loadExternalContent();
    
    // Gestisci i click sul menu laterale
    const sidebarLinks = document.querySelectorAll('.sidebar-menu a');
    
    sidebarLinks.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            
            // Rimuovi classe active da tutti i link
            sidebarLinks.forEach(l => l.classList.remove('active'));
            
            // Aggiungi classe active al link cliccato
            this.classList.add('active');
            
            // Qui puoi caricare contenuti diversi basati sul link cliccato
            const section = this.getAttribute('href').substring(1);
            loadContentForSection(section);
        });
    });
});

// Funzione per caricare contenuti diversi basati sulla sezione
function loadContentForSection(section) {
    const contentContainer = document.getElementById('content-container');
    
    // Mostra loader
    contentContainer.innerHTML = '<div class="loading">Caricamento contenuto...</div>';
    
    // Simula caricamento di contenuti diversi
    setTimeout(() => {
        switch(section) {
            case 'getting-started':
                loadExternalContent();
                break;
            case 'release-notes':
                contentContainer.innerHTML = '<h1>Release Notes</h1><p>Contenuto delle release notes...</p>';
                break;
            case 'tutorials':
                contentContainer.innerHTML = '<h1>Tutorials</h1><p>Contenuto dei tutorial...</p>';
                break;
            default:
                loadExternalContent();
        }
    }, 300);
}

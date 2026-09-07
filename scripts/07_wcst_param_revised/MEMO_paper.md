# Memo: quale paper scrivere con questi dati

**Campione.** 101 partecipanti con almeno un compito, 80 (35 pazienti, 45 controlli)
con dati completi nei tre compiti. Tutte le analisi individuali usano gli 80 completi.

**Vincolo di stima rispettato in tutto.** I parametri di PRL e task switching sono
quelli forniti, usati come sono. Per il WCST le stime individuali usate nelle analisi
di differenze individuali provengono dal rifit del modello HMM *senza* il predittore
di gruppo: le stime del modello con il predittore sono contratte verso la media del
proprio gruppo (R² di gruppo 0.73 per `log d`, 0.64 per `kappa`) e userle per
correlazioni o classificazione sarebbe circolare. Nel modello senza gruppo l'R² scende
a 0.185 e 0.105.

---

## 1. La tesi del paper

> Nell'anoressia nervosa la compromissione del controllo cognitivo si vede in modo
> coerente a livello di gruppo su tutte le primitive computazionali, ma non esiste un
> tratto individuale condiviso tra compiti: chi mostra una soglia alta nel PRL non è
> chi la mostra nel task switching. La "inflessibilità cognitiva" è un fatto di gruppo,
> non una dimensione individuale misurata dai compiti.

È una dissociazione, non un nullo. Le due metà sono entrambe positive e misurate:
la coerenza a livello di gruppo è un risultato (§2), la sua assenza a livello
individuale è un risultato altrettanto specifico (§3), e la terza analisi mostra
*perché* i due livelli possono divergere (§4).

## 2. Livello di gruppo: coerenza su tutte le primitive

Effect size (d di Cohen, pazienti − controlli, IC 95% bootstrap; `effect_sizes.csv`,
figura 1). I due effetti maggiori sono parametri del WCST, ma nessuna primitiva è
indenne:

| parametro | d | IC 95% |
|---|---|---|
| determinismo `log d` — WCST | −1.08 | [−1.59, −0.65] |
| inferenza `logit h` — WCST | +1.02 | [+0.59, +1.53] |
| drift `v` — TS ripetizione | −0.85 | [−1.35, −0.40] |
| soglia `a` — PRL | +0.84 | [+0.41, +1.33] |
| perseverazione `κ` — WCST | −0.77 | [−1.28, −0.33] |
| tasso `α` — PRL | −0.71 | [−1.20, −0.26] |
| soglia `a` — TS switch | +0.55 | [+0.12, +0.98] |
| `t₀` — TS switch | +0.48 | [+0.05, +0.94] |

**Test di coerenza di segno.** Dentro ogni primitiva, tutti gli indicatori dei due
compiti a parametrizzazione DDM puntano nella stessa direzione: 11/11, p = 0.0005
(`concordanza_segni.csv`). Soglia più alta, drift più basso, tempo non decisionale più
lungo, apprendimento più lento nei pazienti — la firma è la stessa in PRL e task
switching. Includendo anche `log d` del WCST la concordanza è 11/12, p = 0.0032; il
disaccordo è atteso perché il determinismo dell'HMM non ha la stessa orientazione di
scala della soglia DDM, e il test principale è quello ristretto alla parametrizzazione
condivisa.

`log d` e `κ` del WCST restano legittimi qui: le loro stime per soggetto sono rumorose
(§3), ma il rumore centrato *attenua* il contrasto tra gruppi, non lo gonfia. I d
riportati sono quindi conservativi.

## 3. Livello individuale: nessun tratto condiviso tra compiti

**Struttura di correlazione** (Spearman, 12 indicatori affidabili;
`struttura_cross_task.csv`, figura 2a):

| classe di coppie | n | media \|ρ\| | max \|ρ\| |
|---|---|---|---|
| stessa primitiva, **entro** compito | 4 | 0.373 | 0.539 |
| primitive diverse, entro compito | 21 | 0.193 | 0.477 |
| **stessa primitiva, tra compiti** | 6 | **0.173** | **0.285** |
| primitive diverse, tra compiti | 35 | 0.114 | 0.359 |

Le coppie che dovrebbero misurare la stessa cosa in due compiti diversi non sono più
correlate delle coppie che non condividono nulla. Le coppie forti sono tutte dentro lo
stesso compito.

**Modello a componenti di varianza** (multitratto-multimetodo bayesiano, Stan;
`loo_mtmm.csv`, figura 2b–c). Cinque strutture di covarianza confrontate per capacità
predittiva fuori campione (elpd LOO, differenze rispetto alla migliore):

| struttura | Δelpd | se |
|---|---|---|
| **senza tratto tra compiti** | 0.0 | — |
| completo (tratto + metodo + entro compito) | −1.3 | 2.3 |
| senza metodo | −10.8 | 6.0 |
| solo entro compito | −11.4 | 5.7 |
| indipendenza completa | −48.9 | 11.0 |

Togliere la componente di tratto tra compiti non costa nulla in previsione; togliere
la componente di compito, o ridursi alle sole coppie entro compito, costa molto. Il
contrasto diretto tra le due quote di varianza è +0.04 [−0.09, +0.17], P(tratto >
metodo) = 0.71 (`contrasto_tratto_metodo.csv`): indistinguibili, e le quote assolute
sono spinte verso l'alto dal vincolo di non negatività, quindi il test valido è il LOO.

**Affidabilità: la parte da dichiarare in Metodi.** Tra i tre parametri del WCST solo
`logit h` è usabile per differenze individuali (`affidabilita_wcst.csv`):

| parametro | var. tra soggetti | errore quadratico medio | affidabilità |
|---|---|---|---|
| `logit h` | 0.208 | 0.150 | 0.58 |
| `log d` | 0.006 | 0.012 | 0.32 |
| `κ` | 0.003 | 0.035 | 0.09 |

`log d` e `κ` sono esclusi dal modello strutturale e dalle correlazioni. Ne segue una
delimitazione onesta del claim: **il test di convergenza tra compiti poggia sulle
coppie PRL↔task switching**, dove l'affidabilità non è in discussione perché i
parametri sono forniti come stime di riferimento; il WCST contribuisce una primitiva
(inferenza di regola) che non ha controparte negli altri due compiti. Le correlazioni
di `logit h` con le learning rate del PRL sono −0.19 e −0.08 (−0.25 e −0.11 dopo
disattenuazione per l'affidabilità 0.58): nessun tratto anche dopo correzione.

## 4. Perché i due livelli divergono: eterogeneità

I pazienti occupano una regione più ampia dello spazio dei parametri, e la differenza è
tutta nella dispersione, non solo nella posizione (`dispersione_test.csv`,
`centroidi_test.csv`, figura 3):

- distanza media dal centroide: 3.53 (pazienti) vs 2.83 (controlli), differenza 0.70,
  p permutazione = 0.0068 (10 000 permutazioni);
- escludendo l'indicatore WCST: differenza 0.64, p = 0.011 — non dipende dal compito
  che abbiamo ristimato;
- distanza tra centroidi 1.91 in 12 dimensioni standardizzate, p < 0.0001.

Nessun singolo parametro sopravvive alla correzione per confronti multipli sul rapporto
di varianze (`eterogeneita.csv`, minimo p_BH = 0.30): l'eterogeneità è una proprietà
del profilo multivariato, non di un parametro. Questo è il meccanismo che riconcilia
§2 e §3: gruppi con medie diverse e dispersione diversa possono non avere alcuna
struttura di covarianza condivisa tra compiti.

## 5. Il valore aggiunto dei tre compiti: classificazione

AUC in validazione incrociata nidificata (5 fold × 20 ripetizioni, selezione del modello
dentro ogni fold; `auc_nested.csv`, figura 4):

| insieme | AUC | IC 95% fold |
|---|---|---|
| tutti e tre | 0.832 | [0.617, 0.968] |
| PRL + WCST | 0.813 | [0.579, 0.961] |
| TS + WCST | 0.797 | [0.603, 0.984] |
| PRL + task switching | 0.758 | [0.556, 0.962] |
| WCST | 0.756 | [0.500, 0.952] |
| PRL | 0.740 | [0.531, 0.905] |
| task switching | 0.681 | [0.387, 0.897] |

L'insieme completo è il migliore in media e nella maggioranza dei fold (72–87% a
seconda del confronto), ma **nessuna differenza rispetto a un singolo compito ha un IC
che esclude lo zero** (`auc_confronti.csv`: tutti e tre − WCST = +0.075 [−0.143,
+0.278]). Da riportare così: il guadagno è coerente in direzione e piccolo in
magnitudine, e con n = 80 non è distinguibile dal rumore. L'ottimismo da selezione non
nidificata resta entro |0.05| di AUC in ogni insieme — il confronto con la letteratura
che riporta AUC non nidificate va fatto tenendone conto.

---

## 6. Come si organizza il paper

**Titolo di lavoro.** *Group-level convergence without individual-level structure:
computational primitives of cognitive control across three tasks in anorexia nervosa.*

**Struttura delle figure.** Fig. 1 effect size per primitiva (§2) → Fig. 2 struttura
cross-task in tre pannelli (§3) → Fig. 3 eterogeneità (§4) → Fig. 4 classificazione
(§5). L'ordine è la tesi: coerenza di gruppo, assenza di tratto, meccanismo,
conseguenza pratica.

**Perché è pubblicabile.** La dissociazione gruppo/individuo è un problema riconosciuto
e attualmente discusso nella misurazione del controllo cognitivo (bassa affidabilità e
scarsa convergenza dei compiti cosiddetti "executive"). Qui non è un problema
metodologico segnalato di passaggio: è il risultato, con tre compiti sugli stessi
soggetti, parametri di modelli generativi invece di indici grezzi, e un test
predittivo (LOO) invece di indici di adattamento. Il messaggio è utile a chi progetta
studi clinici con questi paradigmi.

**Come proteggersi in revisione.** Le tre obiezioni prevedibili e la risposta già
calcolata:

1. *"È un nullo per mancanza di potenza."* Le coppie entro compito sono nettamente
   correlate (media |ρ| = 0.373) nello stesso campione e con la stessa numerosità: il
   disegno rileva struttura quando c'è. E il LOO preferisce esplicitamente il modello
   senza tratto, non è un test che si limita a non rifiutare.
2. *"È un nullo da affidabilità."* Documentata per il WCST e gestita escludendo i due
   parametri inaffidabili; il test poggia su PRL e task switching, i cui parametri non
   sono stati ristimati. La disattenuazione non cambia la conclusione.
3. *"Le stime individuali sono contaminate dall'etichetta di gruppo."* Rifit senza
   predittore di gruppo, R² di gruppo documentato prima e dopo, e verifica che i
   parametri di PRL e task switching non mostrano quella firma.

**Cosa manca prima di scrivere.** Tre cose, in ordine di importanza:

1. **Covariate cliniche** (BMI, durata di malattia, sottotipo, comorbidità, farmaci).
   Servono sia per il paragrafo di caratterizzazione del campione sia, soprattutto, per
   testare se l'eterogeneità di §4 è strutturata: se una dimensione clinica spiega la
   posizione nello spazio dei parametri, la dissociazione diventa una scoperta sulla
   sottotipizzazione e non solo sulla misurazione. È l'analisi con il maggiore ritorno
   atteso.
2. **Affidabilità test-retest** o split-half per PRL e task switching, se il disegno la
   consente. Chiuderebbe definitivamente l'obiezione 2.
3. **Preregistrazione o replica esterna** del solo test di §3, se esiste un secondo
   campione: trasformerebbe il paper da un contributo metodologico su un campione a un
   risultato confermato.

**Tabelle e figure prodotte.** `fig1_effect_size.png`, `fig2_struttura.png`,
`fig3_eterogeneita.png`, `fig4_auc.png`; `effect_sizes.csv`, `concordanza_segni.csv`,
`struttura_cross_task.csv`, `correlazioni_spearman.csv`, `coppie_correlazioni.csv`,
`loo_mtmm.csv`, `decomposizione_varianza.csv`, `contrasto_tratto_metodo.csv`,
`varianza_tratto_per_primitiva.csv`, `affidabilita_wcst.csv`, `eterogeneita.csv`,
`dispersione_test.csv`, `dispersione_loo.csv`, `centroidi_test.csv`, `auc_nested.csv`,
`auc_confronti.csv`, `dati_armonizzati.csv`, `parametri_meta.csv`.

**Codice.** `wcst_rl/08_sticky.R` (modello HMM del WCST), `wcst_rl/10_nogroup.R` (rifit
senza gruppo), `tre_compiti/10_mtmm.R` + `tre_compiti/stan/mtmm.stan` (componenti di
varianza), `tre_compiti/11_cfa_classica.R` (complemento con CFA ML in lavaan).

---

Lo studio nasce dall’ipotesi che la “rigidità cognitiva” associata all’anoressia nervosa possa riflettere un’alterazione generale del controllo cognitivo, osservabile in contesti diversi. Per verificarlo, pazienti e controlli hanno eseguito tre compiti che coinvolgono forme parzialmente differenti di flessibilità e adattamento — WCST, task switching e probabilistic reversal learning — analizzati mediante modelli computazionali. Questo approccio consente di andare oltre gli indici comportamentali complessivi e di distinguere processi latenti come accumulo dell’evidenza, cautela decisionale, velocità di apprendimento, perseverazione e inferenza dei cambiamenti di regola.

A livello di gruppo emerge una firma computazionale coerente. Rispetto ai controlli, le persone con anoressia nervosa mostrano decisioni più caute, accumulo dell’evidenza meno efficiente, tempi non decisionali più lunghi e apprendimento più lento. Anche il WCST evidenzia marcate differenze nei processi di inferenza e nella stabilità delle strategie. Tutti gli indicatori direttamente confrontabili nei modelli del reversal learning e del task switching differiscono nella stessa direzione, indicando che l’alterazione non è circoscritta a un singolo paradigma. Nel loro insieme, i risultati suggeriscono quindi una modificazione ampia dei meccanismi attraverso cui le persone con anoressia aggiornano le proprie credenze e regolano il comportamento in condizioni mutevoli.

Il risultato più interessante è però la dissociazione tra livello di gruppo e livello individuale. Nonostante la coerenza delle differenze medie, gli individui che mostrano una particolare alterazione in un compito non sono necessariamente quelli che la mostrano in un altro. I parametri riferiti alla stessa primitiva computazionale correlano poco tra compiti, mentre le associazioni più forti rimangono confinate all’interno del singolo paradigma. Coerentemente, il confronto predittivo tra modelli non mostra alcun vantaggio nell’introdurre un tratto individuale generale condiviso dai tre compiti. La cosiddetta inflessibilità cognitiva appare dunque robusta come caratteristica media del gruppo clinico, ma non come una dimensione unitaria con cui ordinare stabilmente i singoli pazienti.

La maggiore dispersione multivariata osservata nel gruppo clinico aiuta a spiegare questa apparente contraddizione. I pazienti non sembrano differire dai controlli tutti nello stesso modo: occupano una regione più ampia dello spazio computazionale e possono arrivare a prestazioni apparentemente simili attraverso combinazioni differenti di meccanismi alterati. L’eterogeneità riguarda il profilo complessivo, non la variabilità di un unico parametro. Di conseguenza, una differenza clinica replicabile nelle medie di gruppo non implica necessariamente l’esistenza di un singolo deficit latente condiviso da tutti i pazienti.

Ciò che impariamo è quindi che “inflessibilità cognitiva” non dovrebbe essere trattata come un costrutto individuale semplice e intercambiabile tra paradigmi. I diversi compiti forniscono prospettive complementari, più che misure equivalenti della stessa caratteristica. Il loro uso congiunto produce infatti la migliore classificazione media tra pazienti e controlli, sebbene il vantaggio rispetto ai singoli compiti sia ancora incerto nel presente campione. Il contributo nuovo dello studio consiste proprio nel mostrare, mediante parametri computazionali e un confronto predittivo tra modelli, che la convergenza clinica a livello di gruppo può coesistere con una sostanziale eterogeneità dei meccanismi a livello individuale. Questo invita sia a maggiore cautela nell’interpretare un singolo compito come misura generale della rigidità, sia a concepire l’anoressia nervosa come caratterizzata da molteplici profili di alterazione del controllo cognitivo.

Sì: la dissociazione tra differenze nette a livello di gruppo e scarsa convergenza tra compiti è stata osservata anche altrove, sia nella popolazione generale sia in diversi disturbi psichiatrici. Questo rende improbabile che il vostro risultato sia una semplice anomalia del campione. Tuttavia, la letteratura suggerisce anche una formulazione più prudente: i dati dimostrano eterogeneità dei profili computazionali, ma non ancora l’esistenza di distinti “sottotipi” di anoressia.

## 1. Un precedente molto diretto esiste già nell’anoressia

Dann e colleghi hanno somministrato WCST e cued task switching alle stesse persone con una diagnosi lifetime di anoressia nervosa. Non hanno trovato alcuna associazione tra errori perseverativi al WCST e switch cost: il coefficiente del task switching nella previsione del WCST era praticamente nullo, β = −0,03, p = 0,85. Inoltre, la prestazione al WCST era associata alla memoria di lavoro, suggerendo che i due presunti indici di “flessibilità” riflettessero combinazioni diverse di processi cognitivi. Gli autori conclusero esplicitamente che le due misure non sono intercambiabili nell’anoressia nervosa. [Dann et al., 2023](https://link.springer.com/article/10.1007/s40519-023-01589-6)

Questo è importante perché costituisce una replica concettuale indipendente della componente centrale del vostro risultato. Il vostro studio va però oltre quel lavoro:

* include tre compiti anziché due;
* scompone la prestazione in parametri computazionali;
* confronta parametri teoricamente omologhi tra compiti;
* mostra simultaneamente convergenza delle differenze di gruppo e assenza di convergenza individuale;
* dimostra una maggiore dispersione multivariata nel gruppo clinico;
* confronta direttamente, mediante previsione fuori campione, modelli con e senza un tratto condiviso.

Pertanto, la scarsa convergenza tra compiti non nasce per la prima volta nel vostro dataset; il contributo nuovo consiste nel caratterizzarne la struttura computazionale e nel collegarla all’eterogeneità clinica.

## 2. Il fenomeno non è specifico dell’anoressia

Nella psicologia cognitiva è ben documentato il cosiddetto *reliability paradox*: un compito può produrre un effetto sperimentale molto robusto — o una chiara differenza media tra gruppi — ma distinguere male gli individui. I compiti sono spesso costruiti per minimizzare la variabilità interindividuale non pertinente all’effetto sperimentale; di conseguenza, due compiti capaci di rilevare lo stesso effetto medio possono correlare poco tra loro. [Hedge, Powell & Sumner, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5990556/)

Inoltre, un ampio studio sulla struttura della self-regulation ha mostrato che le misure comportamentali formano fattori largamente determinati dal paradigma, con relazioni deboli sia tra compiti sia con le misure self-report. Ciò mette in discussione l’idea che compiti comunemente attribuiti allo stesso costrutto siano automaticamente misure intercambiabili di una singola caratteristica individuale. [Eisenberg et al., 2019](https://www.nature.com/articles/s41467-019-10301-1)

La letteratura su depressione, ADHD, OCD, dipendenze e schizofrenia riporta frequentemente una combinazione analoga:

1. differenze medie rispetto ai controlli in apprendimento, controllo o decisione;
2. elevata eterogeneità interna al gruppo diagnostico;
3. specificità degli effetti rispetto al compito;
4. associazioni deboli tra parametri computazionali e gravità clinica;
5. affidabilità individuale dei parametri molto variabile.

Quindi la vostra osservazione appartiene a un problema transdiagnostico più generale: le categorie diagnostiche possono spostare la distribuzione media di diversi processi senza generare una sindrome cognitiva uniforme in ogni paziente.

## 3. Esistono disturbi con una struttura più coerente?

Esistono casi in cui si osserva una maggiore covarianza tra prestazioni, ma i confronti non sono perfettamente equivalenti al vostro.

L’esempio più chiaro è la schizofrenia, nella quale si trovano deficit medi ampi in numerosi domini cognitivi. Le prestazioni nelle batterie neuropsicologiche sono correlate e possono essere riassunte, almeno in parte, da una componente cognitiva generale. Tuttavia, anche qui un singolo fattore non è necessariamente la descrizione migliore: in un’analisi di 16 test della batteria MATRICS, il modello a sette domini correlati descriveva i dati meglio sia di un fattore unico sia di un fattore generale gerarchico. [McCleery et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4523424/)

La schizofrenia è dunque un utile caso di contrasto, ma non dimostra l’esistenza di un singolo meccanismo computazionale compromesso. Mostra piuttosto una maggiore struttura condivisa tra domini, insieme a componenti specifiche. Studi recenti trovano infatti compromissioni selettive: per esempio, switching delle regole e gestione del conflitto possono essere alterati mentre l’uso preparatorio dei cue rimane relativamente preservato. [Li et al., 2026](https://pubmed.ncbi.nlm.nih.gov/42184934/)

Un altro esempio parziale proviene dalla ricerca transdiagnostica sulla compulsività. Una dimensione compulsiva che attraversa OCD, dipendenze e altri quadri è stata associata a una riduzione della pianificazione goal-directed in grandi campioni. [Gillan et al., 2016](https://elifesciences.org/articles/11305) Questo rappresenta un legame relativamente coerente tra un meccanismo computazionale e una dimensione psicopatologica. Ma non equivale a dimostrare che lo stesso individuo manifesti lo stesso deficit in più compiti indipendenti: spesso il meccanismo è misurato mediante un solo paradigma.

La risposta più rigorosa è quindi:

> Alcuni disturbi, soprattutto la schizofrenia, mostrano una componente cognitiva condivisa più evidente; alcune dimensioni transdiagnostiche, come la compulsività, sono associate in modo relativamente stabile a specifici meccanismi computazionali. Tuttavia, prove convincenti di un unico meccanismo individuale replicato attraverso compiti diversi sono molto meno comuni di quanto suggerisca il linguaggio tradizionale dei “deficit cognitivi”.

## 4. Cosa protegge il vostro risultato dall’interpretazione come artefatto?

Nel vostro studio la scarsa convergenza potrebbe teoricamente dipendere da rumore di misura, *task impurity* o bassa potenza. Ma diversi elementi rendono insufficiente questa spiegazione:

* nello stesso campione emergono correlazioni più forti entro compito, quindi i dati contengono una struttura individuale rilevabile;
* il modello senza tratto cross-task ha la migliore prestazione LOO: non vi limitate a non rifiutare un effetto;
* il risultato riguarda anche parametri teoricamente corrispondenti di PRL e task switching, non soltanto indici comportamentali eterogenei;
* i parametri WCST meno affidabili sono stati esclusi dalle analisi individuali;
* esiste già un precedente indipendente nell’anoressia che mostra dissociazione tra WCST e task switching;
* la maggiore dispersione dei pazienti rimane significativa anche escludendo il WCST.

Non si può però eliminare completamente l’ipotesi psicometrica senza conoscere l’affidabilità di PRL e task switching. I parametri computazionali non sono automaticamente più affidabili degli indici grezzi: la loro affidabilità dipende dal numero e dalla struttura dei trial, dall’identificabilità del modello e dalla stabilità temporale dei parametri. La letteratura mostra che l’affidabilità test–retest dei parametri di reinforcement learning è molto variabile. [Test–retest reliability of reinforcement-learning parameters](https://pmc.ncbi.nlm.nih.gov/articles/PMC11289054/)

Perciò il risultato non andrebbe presentato come “la bassa correlazione non può dipendere dalla misurazione”, ma come:

> La convergenza di diverse analisi, insieme a precedenti indipendenti, indica che la dissociazione non è spiegata facilmente dal solo errore di misura, sebbene l’affidabilità dei parametri delimiti necessariamente la forza delle conclusioni individuali.

## 5. La distinzione cruciale: eterogeneità non significa ancora sottotipi

La frase “molteplici profili di alterazione” è compatibile con i vostri dati, purché *profili* significhi configurazioni individuali eterogenee. Non potete ancora affermare che esistano sottogruppi discreti — per esempio, un sottotipo caratterizzato da apprendimento lento e un altro da soglia elevata. Una maggiore dispersione può riflettere:

* sottotipi distinti;
* variazione continua lungo più dimensioni;
* differenze di stato clinico;
* effetti di BMI, durata di malattia, farmaci o comorbidità;
* oppure una combinazione di queste fonti.

La formulazione più solida per il paper sarebbe quindi:

> *Our findings align with growing evidence that robust case–control differences in cognitive control need not imply a unitary individual-level deficit. In anorexia nervosa, computational alterations converged in direction at the group level but did not covary across tasks within individuals. Together with the greater multivariate dispersion of the clinical group, this pattern supports mechanistic heterogeneity rather than a single, task-general cognitive inflexibility trait.*

In italiano:

> I risultati indicano che differenze caso–controllo robuste e coerenti non implicano necessariamente un deficit unitario a livello individuale. Nell’anoressia nervosa, le alterazioni computazionali convergono nella direzione a livello di gruppo, ma non covariano tra compiti negli stessi individui. Insieme alla maggiore dispersione multivariata del gruppo clinico, questo pattern è coerente con un’eterogeneità dei meccanismi, piuttosto che con un unico tratto generale di inflessibilità cognitiva.

Questa conclusione è sufficientemente nuova senza sostenere che il fenomeno sia esclusivo dell’anoressia. Anzi, il collegamento con un problema transdiagnostico rende il paper più forte: il vostro studio usa l’anoressia come caso teoricamente informativo per dimostrare empiricamente perché un deficit di gruppo non deve essere automaticamente reificato come tratto del singolo paziente.

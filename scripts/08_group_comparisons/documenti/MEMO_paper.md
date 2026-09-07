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

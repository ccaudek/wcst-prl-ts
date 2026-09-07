# WCST: perché i parametri non differenziavano i soggetti, e il modello che li differenzia

Tutti i numeri di questo documento sono ricalcolati dai dati grezzi che hai allegato.
La pipeline che li produce è in `wcst_rl/` ed è eseguibile dall'inizio alla fine.

---

## 1. La causa non era lo shrinkage né il numero di trial

La variabile di scelta usata per stimare il modello non conteneva le scelte dei
partecipanti. Nella pipeline originale (`funs_input_for_stan_wcst.R`) il campo
`resp_choice` viene costruito dalla carta *mostrata*, non dalla pila *scelta*:
la matrice `resp_choice` in `wcst_stan_data.RDS` è **identica per tutti gli 88
soggetti** in tutte le 60 colonne. Il modello riceveva quindi la stessa
sequenza di risposte da ogni partecipante, e le sole differenze individuali che
poteva vedere erano quelle nel feedback.

Che non fosse un problema di stima si vede dal confronto sui dati corretti: con
la variabile giusta (`chosen_card`) le scelte coincidono con quelle vecchie solo
nel 25% dei trial, e la coerenza scelta/feedback ricostruita è 1.000 (prima era
incoerente). Vedi `fig_diagnosi_wcst.png`.

Il numero di trial (60) e la gerarchia non erano il problema: erano il rimedio.

## 2. Il modello finale

`wcst_rl/stan/hmm_sticky.stan` — inferenza bayesiana della regola con
perseverazione sulla dimensione:

| componente | parametro | ruolo |
|---|---|---|
| tasso di cambio della credenza | `h` | quanto rapidamente il soggetto abbandona l'ipotesi corrente; per soggetto |
| consistenza decisionale | `d` | quanto la scelta segue la credenza; per soggetto |
| perseverazione | `kappa` | bonus alla dimensione seguita al trial precedente; per soggetto |
| lapse | `lapse` | risposte casuali; solo a livello di popolazione |
| rumore sul feedback | `eta` | probabilità di leggere male il feedback; solo a livello di popolazione |

Struttura gerarchica non centrata, effetto di gruppo stimato *dentro* il
modello (`bgrp`), prior LKJ(2) sulla correlazione tra i tre parametri
individuali. `lapse` ed `eta` restano a livello di popolazione perché il
modello a quattro parametri variabili (`hmm_hier4.stan`) non guadagna nulla
(elpd -1.8 ± 0.8) e lascia le loro SD tra soggetti determinate dal prior.

Valori di popolazione: h = 0.149, d = 5.55, kappa = 0.95, lapse = 0.002,
eta = 0.002.

**Un vincolo necessario.** `eta` va vincolato sotto 0.5: sopra quella soglia il
feedback viene letto invertito e il modello ha un secondo modo esattamente
equivalente. Senza il vincolo, R-hat = 1.5 ed ESS = 7. Con il vincolo,
R-hat massimo 1.007, 0 divergenze.

## 3. Il modello è adeguato

Confronto LOO (elpd per trial; il caso vale -1.386):

| modello | elpd/trial | differenza vs finale | Pareto k > 0.7 |
|---|---|---|---|
| RW sulle pile (pipeline originale) | ~ caso | — | — |
| RW sulle dimensioni + perseverazione | -0.444 | -679 ± 43 | 90 |
| HMM, 4 parametri variabili | -0.316 | -65.6 | 0 |
| HMM, 2 parametri variabili | -0.316 | -63.8 ± 12.2 | 0 |
| **HMM + perseverazione** | **-0.304** | — | **0** |

Il RW sulle *pile* non può superare il caso per un motivo strutturale, non
statistico: l'identità della pila non porta informazione sulla regola, perché
la corrispondenza carta-pila cambia a ogni trial. Era il modello della pipeline
originale.

Posterior predictive check su 8 firme comportamentali classiche, simulando in
avanti contro lo schedule reale del compito: 5 su 8 dentro l'intervallo
predittivo al 90%. Le tre che restano fuori sono piccole e nella stessa
direzione (i partecipanti sono un po' più "appiccicosi" del modello):
errori non perseverativi 0.064 vs 0.054, win-stay 0.978 vs 0.970,
recupero post-switch 0.547 vs 0.569. A livello individuale le correlazioni
osservato-predetto vanno da 0.50 (recupero post-switch) a 0.93 (accuratezza).

## 4. Quali parametri differenziano i soggetti — e quali no

Questa è la parte da riportare con onestà, perché la risposta è asimmetrica.

| parametro | affidabilità | recovery al netto del gruppo | varianza spiegata dal gruppo |
|---|---|---|---|
| `h` | 0.67 | **0.79** | 23% |
| `log d` | 0.51 | 0.69 | 71% |
| `kappa` | 0.46 | **-0.22** | 63% |

- **`h` è utilizzabile come misura individuale.** Affidabilità modello-based
  0.67; split-half su blocchi alternati (30 trial per metà) r = 0.55, corretta
  con Spearman-Brown 0.71, e 0.41 anche al netto del gruppo. Le stime per
  soggetto vanno da 0.11 a 0.45.
- **`d` è quasi solo un contrasto di gruppo.** Il 71% della sua varianza tra
  soggetti è spiegata dall'appartenenza al gruppo; lo split-half grezzo è 0.93
  ma scende a 0.34 al netto del gruppo. La correlazione grezza alta è un
  artefatto dello shrinkage verso due medie di gruppo, non differenziazione
  individuale.
- **`kappa` serve al modello ma non è misurabile per soggetto.** Migliora la
  predizione (+63.8 in elpd) e assorbe la perseverazione che altrimenti
  contaminerebbe `h`, ma il parameter recovery individuale al netto del gruppo
  è nullo (-0.22). Va tenuto nel modello e non interpretato come indice
  individuale.

Conclusione operativa: con 60 trial il WCST supporta **un** parametro
individuale, non tre. Dichiararlo è più solido che presentare tre punteggi di
cui due non reggono.

## 5. Validità di `h`

Convergente con gli indici classici (Spearman, incertezza propagata draw per draw):
accuratezza -0.68 [-0.77, -0.57], errori non perseverativi +0.58 [0.47, 0.68],
ripetizione della dimensione -0.59 [-0.70, -0.48], errori perseverativi
+0.46 [0.36, 0.56], lose-shift -0.38, recupero post-switch -0.33.

Nota interpretativa importante: `h` alto correla **positivamente** con gli
errori perseverativi classici. È il motivo per cui l'indice classico è
ambiguo — conta come "perseverativi" anche errori che nascono da una credenza
instabile, non da rigidità. Il modello separa le due cose; l'indice no.

Tra compiti: `h` correla -0.31 [-0.42, -0.19] con il drift rate `v_0` del
task switching (n = 83) e -0.24 [-0.36, -0.13] con `v_1`; -0.17 [-0.30, -0.04]
con `alpha` del PRL (n = 81). Legami modesti ma con intervalli che escludono
lo zero: `h` non è un artefatto del singolo compito.

## 6. Il risultato di gruppo, e perché è interessante

| parametro | pazienti AN − controlli | IC 90% | P(direzione) |
|---|---|---|---|
| `logit h` | **+0.40** | [+0.10, +0.72] | 0.98 |
| `log d` | -0.14 | [-0.26, -0.03] | 0.98 |
| `kappa` | -0.12 | [-0.32, +0.10] | 0.81 |

Le pazienti hanno un tasso di cambio della credenza **più alto** (h = 0.224 vs
0.159) e una consistenza decisionale più bassa; **non** hanno più
perseverazione. Il segno conta: la differenza non è nella rigidità, ma nella
stabilità della rappresentazione della regola. Il quadro classico
"AN = perseverazione al WCST" viene riprodotto negli indici comportamentali,
ma il modello mostra che quegli errori perseverativi sono generati da una
credenza che si destabilizza troppo facilmente, non da un attaccamento
eccessivo all'ipotesi corrente.

È questo l'argomento del paper: **il modello riclassifica il fenotipo, non solo
lo quantifica**. E arriva con la contropartita metodologica necessaria — un
solo parametro è misurabile a livello individuale, e lo diciamo.

## 7. Come si esegue la pipeline

```
wcst_rl/01_build_stan_data.R      ricostruzione delle scelte (la correzione)
wcst_rl/02_fit_models.R           stima di hmm_hier, hmm_hier4, rw_dim_hier
wcst_rl/03_loo_compare.R          confronto LOO a livello di trial e di soggetto
wcst_rl/04_reliability.R          affidabilità modello-based + split-half
wcst_rl/05_recovery.R             parameter recovery
wcst_rl/06_ppc.R                  posterior predictive check
wcst_rl/07_external.R             validità esterna (hmm_hier)
wcst_rl/08_sticky.R               modello finale: stima, LOO, PPC, affidabilità, recovery
wcst_rl/09_external_sticky.R      validità esterna e tabelle finali (hmm_sticky)
wcst_rl/funs_simulate.R           simulazione generativa e firme comportamentali
```

Tempi indicativi su 10 core: ogni modello 12-17 minuti (4 catene, 1000+1000),
il PPC ~3 minuti.

## 8. Cose da sistemare prima di sottomettere

1. **Verificare la ricostruzione delle scelte con l'output grezzo di
   PsyToolkit.** Ho ricostruito `chosen_card` dalle colonne del file raw e la
   coerenza scelta/feedback è 1.000, ma la conferma va fatta su un paio di
   soggetti a mano. Se la pipeline originale ha alimentato anche altri lavori,
   il bug li riguarda.
2. **Il LOO per soggetto ha 48/88 Pareto k > 0.7.** Va riportato come
   indicativo; il confronto solido è quello a livello di trial (0 k > 0.7).
3. **`b/sigma` per `log d` vale -1.64**: l'effetto di gruppo è più grande della
   SD residua tra soggetti. Coerente con il punto 4, ma da non presentare come
   "effect size individuale".
4. **Le tre firme fuori intervallo** suggeriscono che manca un pezzo di
   inerzia della risposta (stickiness sulla pila, non sulla dimensione).
   Testabile con una variante; non cambia le conclusioni sui gruppi.
5. **PRL e task switching** vanno rifatti con la stessa verifica: se il bug
   sulle scelte è nella funzione condivisa di lettura dei raw, potrebbe
   toccare anche loro. Da controllare prima di usare i loro parametri in un
   modello congiunto.

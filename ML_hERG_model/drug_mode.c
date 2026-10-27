/*
 * drug_mode.c -- drug injection mode for the hardware cardiac cell.
 *
 * Integrates with the existing reduced Priebe-Beuckelmann simulation.
 * The ONLY thing this does to the model is scale three conductances:
 *     gK  (lumped delayed rectifier)  <- hERG block
 *     gNa (fast sodium)               <- INa  block
 *     gCa (L-type calcium)            <- ICaL block
 * Everything else in the ODE is untouched.
 */

#include <math.h>
#include <stdint.h>
#include "drugs.h"

/* Nominal (drug-free) conductances for the epicardial reduced model.
 * These must match the values your existing firmware uses. */
#define GK_NOMINAL   0.019f
#define GNA_NOMINAL  16.0f
#define GCA_NOMINAL  0.064f

/* --- UI state ------------------------------------------------------- */
typedef enum {
    MODE_MANUAL = 0,     /* existing: user sets conductances directly */
    MODE_CONDITION,      /* existing: named presets (hyperkalemia, ...) */
    MODE_DRUG            /* new */
} AppMode;

typedef struct {
    AppMode  mode;
    uint8_t  drug_idx;       /* index into DRUGS[] */
    uint8_t  dose_idx;       /* index into DRUG_DOSES[] */
    float    block_pct[DRUG_N_CH];   /* for the on-screen readout */
} DrugState;

static DrugState g_drug = { MODE_MANUAL, 0, 0, {0.0f, 0.0f, 0.0f} };

/* --- Hill equation --------------------------------------------------- */
/*
 * Fraction of a channel still conducting at concentration D:
 *
 *     phi(D) = 1 / (1 + (D / IC50)^h)
 *
 * D == 0          -> phi = 1   (no block)
 * D == IC50       -> phi = 0.5 (half blocked)
 * IC50 very large -> phi ~ 1   ("no data" entries use 1e12)
 */
static float unblocked_fraction(float D, float ic50, float h)
{
    if (D <= 0.0f || ic50 <= 0.0f) {
        return 1.0f;
    }
    /* Guard the pow(): at ratios this small the block is < 1e-6 anyway,
     * and powf() on a tiny argument is both slow and pointless. */
    float ratio = D / ic50;
    if (ratio < 1e-6f) {
        return 1.0f;
    }
    if (ratio > 1e6f) {
        return 0.0f;
    }
    return 1.0f / (1.0f + powf(ratio, h));
}

/* --- Apply the selected drug ----------------------------------------- */
/*
 * Call this once whenever the user changes drug or dose -- NOT every
 * integration step. The conductances stay fixed for the whole run.
 */
void drug_apply(float *gK, float *gNa, float *gCa)
{
    if (g_drug.mode != MODE_DRUG) {
        *gK  = GK_NOMINAL;
        *gNa = GNA_NOMINAL;
        *gCa = GCA_NOMINAL;
        g_drug.block_pct[0] = 0.0f;
        g_drug.block_pct[1] = 0.0f;
        g_drug.block_pct[2] = 0.0f;
        return;
    }

    const Drug *d = &DRUGS[g_drug.drug_idx];
    float conc = DRUG_DOSES[g_drug.dose_idx] * d->cmax;

    float phi_k  = unblocked_fraction(conc, d->ic50[DRUG_CH_HERG], d->hill[DRUG_CH_HERG]);
    float phi_na = unblocked_fraction(conc, d->ic50[DRUG_CH_INA],  d->hill[DRUG_CH_INA]);
    float phi_ca = unblocked_fraction(conc, d->ic50[DRUG_CH_ICAL], d->hill[DRUG_CH_ICAL]);

    *gK  = GK_NOMINAL  * phi_k;
    *gNa = GNA_NOMINAL * phi_na;
    *gCa = GCA_NOMINAL * phi_ca;

    g_drug.block_pct[DRUG_CH_HERG] = 100.0f * (1.0f - phi_k);
    g_drug.block_pct[DRUG_CH_INA]  = 100.0f * (1.0f - phi_na);
    g_drug.block_pct[DRUG_CH_ICAL] = 100.0f * (1.0f - phi_ca);
}

/* --- Live capture check ---------------------------------------------- */
/*
 * At high doses the cell stops producing a real action potential. Any
 * APD measured then is meaningless, so detect it and say so on screen
 * instead of printing a number. Thresholds mirror the Python model.
 */
int drug_beat_captured(float dVdt_max, float v_peak, float v_rest)
{
    return (dVdt_max > 50.0f) && ((v_peak - v_rest) > 80.0f) && (v_peak > 0.0f);
}

/* --- UI helpers ------------------------------------------------------- */
void drug_next(int delta)
{
    int n = (int)g_drug.drug_idx + delta;
    while (n < 0)            n += DRUG_COUNT;
    while (n >= DRUG_COUNT)  n -= DRUG_COUNT;
    g_drug.drug_idx = (uint8_t)n;
}

void dose_next(int delta)
{
    int n = (int)g_drug.dose_idx + delta;
    while (n < 0)              n += DRUG_N_DOSE;
    while (n >= DRUG_N_DOSE)   n -= DRUG_N_DOSE;
    g_drug.dose_idx = (uint8_t)n;
}

void drug_set_mode(AppMode m) { g_drug.mode = m; }

const Drug *drug_current(void)   { return &DRUGS[g_drug.drug_idx]; }
float drug_current_dose(void)    { return DRUG_DOSES[g_drug.dose_idx]; }
float drug_current_conc(void)    { return DRUG_DOSES[g_drug.dose_idx] * DRUGS[g_drug.drug_idx].cmax; }
const float *drug_block_pct(void){ return g_drug.block_pct; }

/* --- Example screen text ---------------------------------------------
 *
 *   Dofetilide          10x
 *   conc 20.0 nM
 *   hERG  79%  CaL  5%  Na  7%
 *   APD90 +163%   [measured]
 *
 * With a predicted drug the last line becomes e.g.
 *   APD90 +12% (4-31%)  [PREDICTED]
 * so the uncertainty is visible rather than implied.
 * ------------------------------------------------------------------- */

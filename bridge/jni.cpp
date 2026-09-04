#include <jni.h>
#include <string>
#include <string.h>

#include "types.h"
#include "parser.h"
#include "molar.h"
#include "stoich.h"
#include "solution.h"
#include "gas.h"
#include "thermo.h"
#include "equilibrium.h"
#include "acid_base.h"
#include "electrochem.h"
#include "colligative.h"
#include "kinetics.h"
#include "nuclear.h"
#include "balancer.h"

// ── Helper ────────────────────────────────────────────

// Konversi ChemResult → JSON string untuk Kotlin
static std::string result_to_json(const ChemResult& r) {
    if (!r.success) {
        std::string s = "{\"success\":false,\"error\":\"";
        s += (r.error ? r.error : "Unknown error");
        s += "\"}";
        return s;
    }
    char buf[128];
    snprintf(buf, sizeof(buf),
        "{\"success\":true,\"value\":%.10g,\"unit\":\"%s\"}",
        r.value, r.unit ? r.unit : "");
    return std::string(buf);
}

static std::string balancer_to_json(const BalancerResult& r) {
    if (!r.success) {
        std::string s = "{\"success\":false,\"error\":\"";
        s += (r.error ? r.error : "Unknown error");
        s += "\"}";
        return s;
    }
    std::string s = "{\"success\":true,\"coefficients\":[";
    for (int i = 0; i < r.size; i++) {
        if (i > 0) s += ",";
        s += std::to_string(r.coefficients[i]);
    }
    s += "]}";
    return s;
}

#define JNI_FN(name) Java_com_yourapp_engine_ChemEngine_##name
#define JSTR(env, s) (env)->NewStringUTF(s.c_str())

extern "C" {

// ── Molar ─────────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(molarMass)(JNIEnv* env, jobject, jstring formula) {
    const char* f = env->GetStringUTFChars(formula, nullptr);
    auto r = result_to_json(molar_mass(f));
    env->ReleaseStringUTFChars(formula, f);
    return JSTR(env, r);
}

JNIEXPORT jstring JNICALL
JNI_FN(masToMol)(JNIEnv* env, jobject, jstring formula, jdouble mass) {
    const char* f = env->GetStringUTFChars(formula, nullptr);
    auto r = result_to_json(mass_to_mol(f, mass));
    env->ReleaseStringUTFChars(formula, f);
    return JSTR(env, r);
}

JNIEXPORT jstring JNICALL
JNI_FN(molToMass)(JNIEnv* env, jobject, jstring formula, jdouble mol) {
    const char* f = env->GetStringUTFChars(formula, nullptr);
    auto r = result_to_json(mol_to_mass(f, mol));
    env->ReleaseStringUTFChars(formula, f);
    return JSTR(env, r);
}

// ── Stoich ────────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(molToParticles)(JNIEnv* env, jobject, jdouble mol) {
    return JSTR(env, result_to_json(mol_to_particles(mol)));
}

JNIEXPORT jstring JNICALL
JNI_FN(particlesToMol)(JNIEnv* env, jobject, jdouble particles) {
    return JSTR(env, result_to_json(particles_to_mol(particles)));
}

JNIEXPORT jstring JNICALL
JNI_FN(molToVolumeStp)(JNIEnv* env, jobject, jdouble mol) {
    return JSTR(env, result_to_json(mol_to_volume_stp(mol)));
}

JNIEXPORT jstring JNICALL
JNI_FN(percentYield)(JNIEnv* env, jobject, jdouble actual, jdouble theoretical) {
    return JSTR(env, result_to_json(percent_yield(actual, theoretical)));
}

// ── Acid Base ─────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(calcPH)(JNIEnv* env, jobject, jdouble hplus) {
    return JSTR(env, result_to_json(calc_pH(hplus)));
}

JNIEXPORT jstring JNICALL
JNI_FN(pHStrongAcid)(JNIEnv* env, jobject, jdouble conc) {
    return JSTR(env, result_to_json(pH_strong_acid(conc)));
}

JNIEXPORT jstring JNICALL
JNI_FN(pHStrongBase)(JNIEnv* env, jobject, jdouble conc) {
    return JSTR(env, result_to_json(pH_strong_base(conc)));
}

JNIEXPORT jstring JNICALL
JNI_FN(pHWeakAcid)(JNIEnv* env, jobject, jdouble Ka, jdouble conc) {
    return JSTR(env, result_to_json(pH_weak_acid(Ka, conc)));
}

JNIEXPORT jstring JNICALL
JNI_FN(pHWeakBase)(JNIEnv* env, jobject, jdouble Kb, jdouble conc) {
    return JSTR(env, result_to_json(pH_weak_base(Kb, conc)));
}

JNIEXPORT jstring JNICALL
JNI_FN(pHBufferAcid)(JNIEnv* env, jobject, jdouble Ka, jdouble molAcid, jdouble molBase) {
    return JSTR(env, result_to_json(pH_buffer_acid(Ka, molAcid, molBase)));
}

// ── Gas ───────────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(idealGasP)(JNIEnv* env, jobject, jdouble n, jdouble T, jdouble V) {
    return JSTR(env, result_to_json(ideal_gas_pressure(n, T, V)));
}

JNIEXPORT jstring JNICALL
JNI_FN(idealGasV)(JNIEnv* env, jobject, jdouble n, jdouble T, jdouble P) {
    return JSTR(env, result_to_json(ideal_gas_volume(n, T, P)));
}

JNIEXPORT jstring JNICALL
JNI_FN(celsiusToKelvin)(JNIEnv* env, jobject, jdouble c) {
    return JSTR(env, result_to_json(celsius_to_kelvin(c)));
}

// ── Thermo ────────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(gibbsEnergy)(JNIEnv* env, jobject, jdouble dH, jdouble dS, jdouble T) {
    return JSTR(env, result_to_json(gibbs_energy(dH, dS, T)));
}

JNIEXPORT jstring JNICALL
JNI_FN(heatTransfer)(JNIEnv* env, jobject, jdouble m, jdouble c, jdouble dT) {
    return JSTR(env, result_to_json(heat_transfer(m, c, dT)));
}

// ── Electrochem ───────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(cellPotential)(JNIEnv* env, jobject, jdouble cathode, jdouble anode) {
    return JSTR(env, result_to_json(cell_potential(cathode, anode)));
}

JNIEXPORT jstring JNICALL
JNI_FN(nernst25)(JNIEnv* env, jobject, jdouble E0, jdouble n, jdouble Q) {
    return JSTR(env, result_to_json(nernst_25(E0, n, Q)));
}

JNIEXPORT jstring JNICALL
JNI_FN(faradayMass)(JNIEnv* env, jobject, jdouble I, jdouble t, jdouble Ar, jint n) {
    return JSTR(env, result_to_json(faraday_mass(I, t, Ar, n)));
}

// ── Colligative ───────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(boilingPointElevation)(JNIEnv* env, jobject, jdouble Kb, jdouble m, jdouble i) {
    return JSTR(env, result_to_json(boiling_point_elevation(Kb, m, i)));
}

JNIEXPORT jstring JNICALL
JNI_FN(freezingPointDepression)(JNIEnv* env, jobject, jdouble Kf, jdouble m, jdouble i) {
    return JSTR(env, result_to_json(freezing_point_depression(Kf, m, i)));
}

JNIEXPORT jstring JNICALL
JNI_FN(osmoticPressure)(JNIEnv* env, jobject, jdouble M, jdouble T, jdouble i) {
    return JSTR(env, result_to_json(osmotic_pressure(M, T, i)));
}

// ── Kinetics ──────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(halfLifeFirst)(JNIEnv* env, jobject, jdouble k) {
    return JSTR(env, result_to_json(half_life_first(k)));
}

JNIEXPORT jstring JNICALL
JNI_FN(concFirstOrder)(JNIEnv* env, jobject, jdouble c0, jdouble k, jdouble t) {
    return JSTR(env, result_to_json(conc_first_order(c0, k, t)));
}

JNIEXPORT jstring JNICALL
JNI_FN(activationEnergy)(JNIEnv* env, jobject,
    jdouble k1, jdouble T1, jdouble k2, jdouble T2) {
    return JSTR(env, result_to_json(activation_energy(k1, T1, k2, T2)));
}

// ── Nuclear ───────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(remainingAmount)(JNIEnv* env, jobject,
    jdouble initial, jdouble halfLife, jdouble time) {
    return JSTR(env, result_to_json(remaining_amount(initial, halfLife, time)));
}

JNIEXPORT jstring JNICALL
JNI_FN(timeElapsed)(JNIEnv* env, jobject,
    jdouble initial, jdouble remaining, jdouble halfLife) {
    return JSTR(env, result_to_json(time_elapsed(initial, remaining, halfLife)));
}

// ── Balancer ──────────────────────────────────────────

JNIEXPORT jstring JNICALL
JNI_FN(balanceEquation)(JNIEnv* env, jobject, jstring equation) {
    const char* eq = env->GetStringUTFChars(equation, nullptr);
    auto r = balancer_to_json(balance_equation(eq));
    env->ReleaseStringUTFChars(equation, eq);
    return JSTR(env, r);
}

} // extern "C"
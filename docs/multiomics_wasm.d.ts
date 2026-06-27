/* tslint:disable */
/* eslint-disable */

/**
 * Run all three modalities and a lightweight integration summary from raw bytes.
 * Returns a JSON string of `{ genomics, transcriptomics, epigenomics, integration }`,
 * or `{"error":"..."}` on failure.
 */
export function analyze_all(vcf: Uint8Array, tsv: Uint8Array, bed: Uint8Array): string;

/**
 * Parse and analyse a BED methylation file from raw bytes.
 * Returns a JSON string of `EpigenomicsSummary`, or `{"error":"..."}` on failure.
 */
export function analyze_epigenomics(data: Uint8Array): string;

/**
 * Parse and analyse a VCF file from raw bytes.
 * Returns a JSON string of `GenomicsSummary`, or `{"error":"..."}` on failure.
 */
export function analyze_genomics(data: Uint8Array): string;

/**
 * Parse and analyse an expression-matrix TSV from raw bytes.
 * Returns a JSON string of `TranscriptomicsSummary`, or `{"error":"..."}` on failure.
 */
export function analyze_transcriptomics(data: Uint8Array): string;

/**
 * Returns the crate version string.
 */
export function get_version(): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly analyze_all: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => void;
    readonly analyze_epigenomics: (a: number, b: number, c: number) => void;
    readonly analyze_genomics: (a: number, b: number, c: number) => void;
    readonly analyze_transcriptomics: (a: number, b: number, c: number) => void;
    readonly get_version: (a: number) => void;
    readonly __wbindgen_export: (a: number) => void;
    readonly __wbindgen_add_to_stack_pointer: (a: number) => number;
    readonly __wbindgen_export2: (a: number, b: number) => number;
    readonly __wbindgen_export3: (a: number, b: number, c: number) => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;

# multiomics-wasm

WebAssembly bindings for the multiomics analysis engine.

## Build

```bash
# Install wasm-pack
cargo install wasm-pack

# Build (outputs to wasm/pkg/)
wasm-pack build wasm --target web
```

## Usage (JavaScript)

```javascript
import init, { analyze_all, get_version } from './pkg/multiomics_wasm.js';

await init();
console.log(get_version());

const vcf  = await fetch('sample.vcf').then(r => r.arrayBuffer());
const tsv  = await fetch('expr.tsv').then(r => r.arrayBuffer());
const bed  = await fetch('meth.bed').then(r => r.arrayBuffer());

const result = analyze_all(
  new Uint8Array(vcf),
  new Uint8Array(tsv),
  new Uint8Array(bed)
);
console.log(JSON.parse(result));
```

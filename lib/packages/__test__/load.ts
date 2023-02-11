import { createReadStream } from 'node:fs';
import * as readline from 'node:readline/promises';
import { generateFields } from './parse-matrix-from-csv';

// old /^(?<real>-?(?:[0-9]+)(?:\.[0-9]+)?)(?<imag>(?:\+|\-)(?:[0-9]+)(?:\.[0-9]+)?)i$/);
const regExpComplexNumber = //  /^(?<real>[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?[0-9]+)?)[+-]{1}(?<imag>[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?[0-9]+)?)i$/;
/^(?<real>-?(?:[0-9]+)(?:\.[0-9]+)?)(?<imag>(?:\+|\-)(?:[0-9]+)(?:\.[0-9]+)?)i$/;

const regExp2 =  /^(?<real>[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?[0-9]+)?)[+-]{1}(?<imag>[-+]?(?:[0-9]*[.])?[0-9]+(?:[eE][-+]?[0-9]+)?)i$/;
export async function loadData(fullPath: string, separator: string | RegExp, hasColumnNames = true, hasRowNames = true, optionallyQuoted = true, fp64 = true) {

    let line = 0;
    let firstDataLineParsed = false;
    const reader = readline.createInterface({
        input: createReadStream(fullPath, { encoding: 'utf8' }),
    });

    const collectedFields: number[] = [];

    for await (const text of reader) {
        line++;
        if (!text) {
            continue; // just skip
        }
        if (line === 351) {
            console.log('debug1');
        }
        const fields = Array.from(generateFields({ text, line }, separator, optionallyQuoted));
        // if the fields are only comments skip
        if (fields.every(f => f.type === 'comment')) {
            continue;
        }
        if (firstDataLineParsed === false && hasColumnNames) { // we expect the first non-comment line to have names
            firstDataLineParsed = true;
            continue; // skip, disregard for now
        }
        for (let i = 0; i < fields.length; i++) {
            if (hasRowNames && i === 0) {
                continue;
            }
            if (fields[i].type === 'comment') {
                continue; // skip comments in the data
            }
            // are these complex values?
            if (fields[i].text) {
                const ℂ: null | { groups?: Record<string, string> } = fields[i].text?.match(regExpComplexNumber);
                const C2: null | { groups?: Record<string, string> } = fields[i].text?.match(regExp2);
                C2;
                if (ℂ?.groups) {
                    const real = parseFloat(ℂ.groups['real']);
                    const imag = parseFloat(ℂ.groups['imag']);
                    collectedFields.push(real, imag);
                    continue;
                }
                /*const ℝ: null | { groups?: Record<string, string> } = fields[i].text?.match(/^(?<real>-?(?:[0-9]+)(?:\.[0-9]+)?)$/);
                if (ℝ?.groups) {
                    const real = parseFloat(ℝ.groups['real']);
                    collectedFields.push(real);
                    continue;
                }*/
                throw new Error(`Error at parsing field: ${i + 1}, at line ${line}, ${fields[i].text}.`);
            }
            throw new Error(`Error parsing file at line ${line}`);
        }
    }
    return fp64 ? new Float64Array(collectedFields) : new Float32Array(collectedFields);
}




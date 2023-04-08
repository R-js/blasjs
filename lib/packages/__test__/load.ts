import { createReadStream } from 'node:fs';
import * as readline from 'node:readline/promises';
import { generateFields } from './parse-matrix-from-csv';


const regExpComplexNumber =
    /^(?<real>-?(?:[0-9]+)(?:\.[0-9]+)?)(?<imag>(?:\+|-)(?:[0-9]+)(?:\.[0-9]+)?)i$/;

export async function loadData(
    fullPath: string,
    separator: string | RegExp,
    hasColumnNames = true,
    hasRowNames = true,
    optionallyQuoted = true,
    fp64 = true
): Promise<Float64Array | Float32Array> {

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
                if (ℂ?.groups) {
                    const real = parseFloat(ℂ.groups['real']);
                    const imag = parseFloat(ℂ.groups['imag']);
                    collectedFields.push(real, imag);
                    continue;
                }
                throw new Error(`Error at parsing field: ${i + 1}, at line ${line}, ${fields[i].text}.`);
            }
            throw new Error(`Error parsing file at line ${line}`);
        }
    }
    return fp64 ? new Float64Array(collectedFields) : new Float32Array(collectedFields);
}




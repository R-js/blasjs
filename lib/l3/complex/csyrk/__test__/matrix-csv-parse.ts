
function* generateLines(csv: string) {
    // parse for line endings
    let lineStart = 0;
    let line = 1;
    let i = 0;
    while (i < csv.length) {
        // if we see a \r\n treat it as one
        if (csv[i] === '\r') {
            const text = csv.slice(lineStart, i);
            (yield { text, line }) as void;
            line++;
            if (csv[i + 1] === '\n') {
                i += 2;
            }
            else {
                i += 1;
            }
            lineStart = i;
            continue;
        }
        if (csv[i] === '\n') {
            const text = csv.slice(lineStart, i);
            i++;
            lineStart = i;
            (yield { text, line }) as void;
            line++
            continue;
        }
        i++;
    }
    if (lineStart < csv.length) {
        // one last return
        const text = csv.slice(lineStart, i);
        yield { text, line };
    }
}

function isSeparator(text: string, separator: string | RegExp): boolean {
    return (typeof separator === 'string') ? text?.[0] === separator : text[0].match(separator) !== null;
}

function absorbSeparator(text: string, i: number, separator: string | RegExp) {
    while (isSeparator(text[i], separator) && i < text.length) {
        i++;
    }
    return i; // absorbed till EOL
}

function isQuoted(text: string, i: number) {
    return text[i] === '"' || text[i] === '\'';
}

// text[i] is the position after the first quote
function absorbBetweenQuotes(text: string, i: number) {
    while (i < text.length
        &&
        (
            isQuoted(text, i) === false
            ||
            (isQuoted(text, i) && text[i - 1] === '\\')
        )
    ) {
        i++;
    }
    return i === text.length ? i : i + 1;
}

function absorbField(text: string, i: number, separator: string | RegExp, optionallyQuoted = true) {
    if (optionallyQuoted && isQuoted(text, i)) {// i am looking for a quote end
        return absorbBetweenQuotes(text, i + 1)
    }
    while (false === isSeparator(text[i], separator) && i < text.length) {
        i++;
    }
    return i;
}

function absorbComment(text: string, i: number) {
    if (text[i] === '#') {
        return text.length;
    }
    return i;
}

function absorbWhiteSpace(text: string, i: number) {
    while ('\t '.includes(text[i]) && i < text.length) {
        i++;
    }
    return i;
}

type FieldType = 'field' | 'quoted-field' | 'comment' | 'column-name';

function* generateFields(oneLine: { text: string, line: number }, separator: string | RegExp, optionallyQuoted = true) {
    let i = 0;
    const { text, line } = oneLine;
    let j = i;
    while (i < text.length) {
        i = absorbWhiteSpace(text, i); // absorb white space outside of fields
        j = absorbComment(text, i);
        if (j > i) {
            (yield { text: text.slice(i, j).trim(), line, s: i, e: j, type: 'comment' as FieldType }) as void;
            i = j;
        }
        j = absorbField(text, i, separator, optionallyQuoted);
        if (j > i) {
            (yield { text: text.slice(i, j).trim(), line, s: i, e: j, type: isQuoted(text, i) ? 'quoted-field' : 'field' as FieldType }) as void;
            i = j;
        }
        i = absorbSeparator(text, i, separator);
    }
}


export const matrixC = `# some comments
"","V1","V2","V3","V4"  # these are the column names
"1",  -0.420971975425178+0.893161777654979i  ,  -0.944630171426201+0.480959398648743i,1.13479800637549+1.04382439362939i,0.201891445890486+0.721753165361698i
"2",  -0.944630171426201+0.480959398648743i  ,0.435076217019291-0.679847411395675i,1.1299678199116+0.42211028141134i,-0.958405822236894+0.06475720755297i
"3",   1.13479800637549+1.04382439362939i    ,0.201891445890486+0.721753165361698i,1.60338539440229-1.33376062747731i,1.08955672997013-0.50574534029778i
"4",   1.1299678199116+0.42211028141134i     ,-0.958405822236894+0.06475720755297i,1.08955672997013-0.50574534029778i,-1.38221172857068+0.47685804375437i
`

export const matrixA = `
[,1]                [,2]
[1,], -0.7939294-1.4826152i,  1.965057-0.312665i
[2,],  0.0044153+1.0705172i,  1.719627+0.436114i
[3,],  0.5935871+0.0329470i, -1.452431+0.368372i
[4,],  0.5640821+0.9536697i, -0.912274-1.032309i
`;

export default function* parse(text: string, separator: string | RegExp, hasColumnNames = true, hasRowNames = true, optionallyQuoted = true) {

    let firstDataLineParsed = false;
    const lines = generateLines(text);
    //let columnNames: string[] | undefined;
    for (
        let line = lines.next();
        line?.done !== true;
        line = lines.next()
    ) {

        if (line.value?.text === '') {
            continue; // just skip
        }

        if (line.value !== undefined) {
            // this can happen because "generateLines" returns undefined when it runs through all of its "yields"

            // scan fields
            const fields = Array.from(generateFields(line.value, separator, optionallyQuoted));
            if (fields.every(f => f.type === 'comment')) {
                continue; // skip only comment lines
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
                    const ℂ: null | { groups?: Record<string, string> } = fields[i].text?.match(/^(?<real>-?(?:0|1)\.[0-9]+)(?<imag>(?:\+|\-)(?:0|1)\.[0-9]+)i$/);
                    if (ℂ?.groups) {
                        const real = parseFloat(ℂ.groups['real']);
                        const imag = parseFloat(ℂ.groups['imag']);
                        yield* [real, imag];
                        continue;
                    }
                    const ℝ: null | { groups?: Record<string, string> } = fields[i].text?.match(/^(?<real>-?(?:0|1)\.[0-9]+)$/);
                    if (ℝ?.groups){
                        const real = parseFloat(ℝ.groups['real']);
                        yield real;
                        continue;
                    }
                }
            }
        }
    }
}






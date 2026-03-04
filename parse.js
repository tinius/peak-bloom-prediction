import fs from 'fs'
import { parse } from 'csv-parse/sync'
import _ from 'lodash'

const round = (value, exp) => {
	if (typeof exp === 'undefined' || +exp === 0)
		return Math.round(value);

	value = +value;
	exp = +exp;

	if (isNaN(value) || !(typeof exp === 'number' && exp % 1 === 0))
		return NaN;

	value = value.toString().split('e');
	value = Math.round(+(value[0] + 'e' + (value[1] ? (+value[1] + exp) : exp)));

	value = value.toString().split('e');
	return +(value[0] + 'e' + (value[1] ? (+value[1] - exp) : -exp));
}


const data = parse(fs.readFileSync('outputs/latest_summary.csv'), { columns : true, cast: true })

fs.writeFileSync('dist/latest_summary.json', JSON.stringify(data))

const probs = parse(fs.readFileSync('outputs/bloom_probabilities.csv'), { columns : true, cast: true })

fs.writeFileSync('dist/bloom_probabilities.json', JSON.stringify(probs))

const temps = _(parse(fs.readFileSync('outputs/temp_paths_sample50.csv'), { columns : true, cast : true }))
    .groupBy('path_id')
    .mapValues((arr) => {

        return arr.map(row => {
            return {
                date : row.DATE,
                tempF : round(32 + row.tmean_c*1.8, 1)
            }
        }).filter((row, i) => {
            return i <= 73
        })

    })
    .values()
    .valueOf()

fs.writeFileSync('dist/temps.json', JSON.stringify(temps))